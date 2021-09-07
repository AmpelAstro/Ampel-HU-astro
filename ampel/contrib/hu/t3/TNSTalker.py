#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TNSTalker.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 17.11.2018
# Last Modified Date: 04.09.2019
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

import re
from itertools import islice
from typing import Any, Dict, Generator, Iterable, List, Optional, Tuple, TYPE_CHECKING, Union
from ampel.struct.StockAttributes import StockAttributes
from ampel.types import StockId, UBson
from ampel.abstract.AbsT3Unit import AbsT3Unit, T3Send
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.contrib.hu.t3.ampel_tns import (
    TNSClient,
    TNSFILTERID,
    TNS_BASE_URL_REAL,
    TNS_BASE_URL_SANDBOX
)
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id
from ampel.contrib.hu.t3.tns.TNSToken import TNSToken

if TYPE_CHECKING:
    from ampel.content.JournalRecord import JournalRecord


def chunks(l: Iterable, n: int) -> Generator[List, None, None]:
    source = iter(l)
    while True:
        chunk = list(islice(source, n))
        yield chunk
        if len(chunk) < n:
            break


class TNSTalker(AbsT3Unit):
    """
    Get TNS name if existing, and submit selected candidates.
    
    All candidates loaded by T3 will be submitted - it is assumed that *selection* is done
    by an appropriate T2, which also prepares the submit information.
    T2TNSEval is one such example.

    If submit_tns is true, candidates fulfilling the criteria will be sent to the TNS if:
    - They are not known to the TNS OR
    - They are registered by TNS but under a non-ZTF internal name AND resubmit_tns_nonztf set to True OR
    - They are registered by TNS under a ZTF name AND resubmit_tns_ztf is set to True

    if sandbox is set to True it will try to submit candidates to the TNS sandbox, but this API has been unstable
    and might not work properly.
    """

    # TNS config

    # Bot api key frm TNS
    tns_api_key: NamedSecret[dict]
    # Check for TNS for names even if internal name is known
    get_tns_force: bool = False
    # Submit candidates passing criteria (False gives you a 'dry run')
    submit_tns: bool = True
    # Submit all candidates we have a note in the Journal that we submitted this. Overrides the resubmit entries!!
    submit_unless_journal: bool = False
    # Resubmit candidate submitted w/o the same ZTF internal ID
    resubmit_tns_nonztf: bool = True
    # Resubmit candidates even if they have been added with this name before
    resubmit_tns_ztf: bool = False

    # Submit to TNS sandbox only
    sandbox: bool = True
    # weather journal will go to separate collection.
    ext_journal: bool = True

    # AT report config
    base_at_dict: Dict = {
        "reporting_group_id": "82",    # Should be ampel
        "discovery_data_source_id": "48",
        "reporter": "J. Nordin, V. Brinnel, J. van Santen (HU Berlin), A. Gal-Yam, O. Yaron, S. Schulze (Weizmann) on behalf of ZTF",
        "at_type": "1",
    }
    baseremark: str = "See arXiv:1904.05922 for selection criteria."

    slack_token: Optional[NamedSecret] = None
    slack_channel = "#ztf_tns"
    slack_username = "AMPEL"
    # if you have more than this # of reports, send different files
    max_slackmsg_size = 200

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.client = TNSClient(
            TNS_BASE_URL_SANDBOX if self.sandbox else TNS_BASE_URL_REAL,
            self.logger,
            TNSToken(**self.tns_api_key.get()),
        )
        # maintain a second client to check the real TNS if in sandbox mode
        self.reference_client = TNSClient(
            TNS_BASE_URL_REAL,
            self.logger,
            TNSToken(**self.tns_api_key.get()),
        ) if self.sandbox else self.client

    def search_journal_tns(
        self, tran_view: TransientView
    ) -> Tuple[Optional[str], List[str]]:
        """
        Look through the journal for a TNS name.
        Assumes journal entries came from this unit, that the TNS name is saved as "tnsName"
        and internal names as "tnsInternal"
        """
        tns_name, tns_internals = None, []

        def select(entry: "JournalRecord") -> bool:
            return bool(
                (entry["extra"] is not None and ("tnsInternal" in entry["extra"]))
                and entry["unit"]
                and entry["unit"] == self.__class__.__name__
            )

        if jentries := tran_view.get_journal_entries(tier=3, filter_func=select):
            if jentries[-1]["extra"] is not None:
                tns_name = jentries[-1]["extra"].get("tnsName", None)
            tns_internals = [
                entry["extra"].get("tnsInternal", None)
                for entry in jentries
                if entry["extra"] is not None
            ]

        self.logger.info(
            "Journal search",
            extra={
                "tranId": tran_view.id,
                "tnsName": tns_name,
                "tnsInternals": tns_internals,
            },
        )

        return tns_name, tns_internals

    def search_journal_submitted(self, tran_view: TransientView) -> bool:
        """
        Look through the journal for whether this sender submitted this to TNS.
        Assumes journal entries came from this unit, that the TNS name is saved as "tnsName"
        and tnsSender stores the api key used ('tnsSender': self.tns_api_key')
        """

        def select(entry: "JournalRecord") -> bool:
            return bool(
                (
                    entry["extra"] is not None
                    and (entry["extra"].get("tnsSender") == self.tns_api_key.get()["name"])
                    and "tnsSubmitResult" in entry["extra"]
                )
                and entry["unit"]
                and entry["unit"] == self.__class__.__name__
            )

        # Find the latest tns name (skipping previous)
        if tran_view.get_journal_entries(
            tier=3,
            filter_func=select,
        ):
            self.logger.info("TNS submitted", extra={"tnsSender": self.tns_api_key.get()["name"]})
            return True
        else:
            self.logger.info("Not TNS submitted", extra={"tnsSender": self.tns_api_key.get()["name"]})
            return False

    def _query_tns_names(self, tran_view: TransientView, ra: float, dec: float) -> Tuple[Optional[str], List]:
        """
        query the TNS for names and internals at the position
        of the transient.
        """
        # query the TNS for transient at this position. Note that we check the real TNS for names for compatibility...

        tns_name, tns_internal = self.client.getNames(
            ra=ra,
            dec=dec
        )

        # Skip the AT SN prefix if present
        if tns_name is not None:
            tns_name = re.sub("^AT", "", tns_name)
            tns_name = re.sub("^SN", "", tns_name)

        # be nice and then go
        ztf_name = to_ztf_id(tran_view.id)
        self.logger.info(
            "looking for TNS name in the TNS.",
            extra={
                "ZTFname": ztf_name,
                "ra": ra,
                "dec": dec,
                "tnsName": tns_name,
                "tnsInternals": [tns_internal],
            },
        )
        return tns_name, [tns_internal]

    def _find_tns_tran_names(
        self, tran_view: TransientView
    ) -> Tuple[Optional[str], List[str]]:
        """
        search for TNS name in tran_view.tran_names. If found,
        look in the TNS for internal names and return them
        """

        # First, look if we already registered a name
        tns_name, tns_internals = None, []
        names: List[str] = (
            [str(name) for name in (tran_view.stock["name"] or [])]
            if tran_view.stock
            else []
        )
        for tname in names:

            if "TNS" in tname and (not self.get_tns_force):
                self.logger.info(
                    "found TNS name in tran_names.",
                    extra={"TNSname": tname, "TransNames": names},
                )
                # as TNS to give you the internal names.
                # we remove the 'TNS' part of the name, since this has been
                # added by the TNSMatcher T3, plus we skip the prefix
                # We here assume that the AT/SN suffix is cut
                tns_name = tname.replace("TNS", "")
                # Not using sandbox (only checking wrt to full system). 
                tns_internals, runstatus = self.reference_client.getInternalName(tns_name)

        # be nice with the logging
        ztf_name = to_ztf_id(tran_view.id)
        self.logger.info(
            "looked for TNS name in self.tran_names",
            extra={
                "ZTFname": ztf_name,
                "tnsName": tns_name,
                "tnsInternals": tns_internals,
                "TransNames": names,
            },
        )

        return tns_name, tns_internals


    def find_tns_name(
        self, tran_view: TransientView, ra: float, dec: float
    ) -> Tuple[Optional[str], List[str], Optional[JournalAttributes]]:
        """
        extensive search for TNS names in:
        - tran_view.tran_names (if added by TNSMatcher)
        - the journal of tran_view (if added by this T3)
        - the TNS itself (if no name can be found with the above)

        Returns:
        --------
            tns_name, tns_internals, jup: tns_name, tns_internal, and journal update
        """

        ztf_name = to_ztf_id(tran_view.id)
        self.logger.info("looking for TNS name", extra={"ZTFname": ztf_name})

        # first we look in the journal, this is the cheapest option. If we have
        # a valid name from the journal and if you do not want to look again in
        # the TNS, we are fine. NOTE: in this case you don't return a journal update.
        tns_name, tns_internals = self.search_journal_tns(tran_view)
        self.logger.debug("Found tns name in journal: %s" % (tns_name))
        if (not tns_name is None) and (not self.get_tns_force):
            return tns_name, tns_internals, None

        # second option in case there is no TNS name in the journal: go and look in tran_names
        # and if you don't find any, go and ask TNS again.
        tns_name_new, tns_internals_new = self._find_tns_tran_names(tran_view)
        self.logger.debug(
            "Find tns names added to the ampel name list: %s internal %s"
            % (tns_name_new, tns_internals_new)
        )
        if tns_name_new is None:
            tns_name_new, tns_internals_new = self._query_tns_names(tran_view, ra, dec)
            self.logger.debug(
                "Proper check of tns done, found name %s" % (tns_name_new)
            )
        

        # now, it is possible (if you set self.get_tns_force) that the
        # new TNS name is different from the one we had in the journal. We always
        # use the most recent one. In this case we also create a JournalUpdate
        jup = None
        if not tns_name_new is None:

            # what happen if you have a new name that is different from the old one?
            if tns_name is not None and not tns_name == tns_name_new:
                self.logger.info(
                    "Adding new TNS name",
                    extra={"tnsOld": tns_name, "tnsNew": tns_name_new},
                )

                # create content of journal entry. Eventually
                # update the list with the new internal names if any are found
                jcontent = {"tnsName": tns_name_new}
                if tns_internals_new is not None:
                    tns_internals.extend(tns_internals_new)
                    for tns_int in tns_internals_new:
                        jcontent.update({"tnsInternal": tns_int})

                # create a journalUpdate and update the tns_name as well. TODO: check with JNo
                jup = JournalAttributes(extra=jcontent)
                tns_name = tns_name_new

            elif tns_name is None:
                # Set the new name
                self.logger.info(
                    "Adding first TNS name", extra={"tnsNew": tns_name_new}
                )

                # create content of journal entry. Eventually
                # update the list with the new internal names if any are found
                jcontent = {"tnsName": tns_name_new}
                if tns_internals_new is not None:
                    tns_internals.extend(tns_internals_new)
                    for tns_int in tns_internals_new:
                        jcontent.update({"tnsInternal": tns_int})

                # create a journalUpdate and update the tns_name as well. TODO: check with JNo
                jup = JournalAttributes(extra=jcontent)
                tns_name = tns_name_new
                # tns_internals = tns_internals_new

        # bye!
        return tns_name, tns_internals, jup



    def process(self, gen: Generator[TransientView, T3Send, None]) -> None:
        """
        Loop through transients and check for TNS names and/or candidates to submit
        """

        # Reports to be sent, indexed by the transient view IDs (so that we can check in the replies)
        atreports: Dict[StockId, Dict[str, Any]] = {}

        for tran_view in gen:

            ztf_name = to_ztf_id(tran_view.id)

            # Obtain atdict start from T2 result
            t2result = tran_view.get_t2_result(unit_id="T2TNSEval")
            if not isinstance(t2result, dict):
                raise ValueError("Need to have a TNS atdict started from a suitable T2.")
            # Create the submission dictionary - in case the transient is to be submitted
            atdict = dict(t2result['atdict'])
            atdict.update(self.base_at_dict)
            atdict["internal_name"] = ztf_name

            ra, dec = atdict['ra']['value'], atdict['dec']['value']

            self.logger.info(
                "TNS init dict found", extra={"tranId": tran_view.id, "ztfName": ztf_name}
            )

            # Three things we can find out:
            # - Did this AMPEL channel submit the transient (according to Journal)
            # - Did anyone submit a transient with this ZTF name?
            # - Did anyone submit a transient at this position?

            # Simplest case to check. We wish to submit everything not noted as submitted
            if self.submit_unless_journal:
                if self.search_journal_submitted(tran_view):
                    # Note already submitted
                    self.logger.info("ztf submitted", extra={"ztfSubmitted": True})
                else:
                    # add AT report
                    self.logger.info("Add TNS report list", extra={"id": tran_view.id})
                    atreports[tran_view.id] = atdict
                continue


            # find the TNS name, either from the journal, from tran_names, or
            # from TNS itself. If new names are found, create a new JournalUpdate
            tns_name, tns_internals, jup = self.find_tns_name(tran_view, ra, dec)
            self.logger.debug("TNS got %s internals %s" % (tns_name, tns_internals))

            if tns_name is not None:

                # Chech whether this ID has been submitted (note that we do not check
                # whether the same candidate was submitted as different ZTF name) and
                # depending on what's already on the TNS we can chose to submit or not
                is_ztfsubmitted = ztf_name in tns_internals
                if is_ztfsubmitted:
                    # Already registered under this name. Only submit if we explicitly configured to do this
                    if not self.resubmit_tns_ztf:
                        self.logger.info(
                            "ztf submitted",
                            extra={
                                "ztfSubmitted": is_ztfsubmitted,
                                "tnsInternals": tns_internals,
                            },
                        )
                        continue

                # Also allow for the option to not submit if someone (anyone) already did this. Not sure why this would be a good idea.
                if not is_ztfsubmitted and not self.resubmit_tns_nonztf:
                    self.logger.info(
                        "already in tns, skipping",
                        extra={
                            "ztfSubmitted": is_ztfsubmitted,
                            "tnsInternals": tns_internals,
                        },
                    )
                    continue

            # Passed all cuts, add to submit list
            self.logger.info("Added to report list")
            atreports[tran_view.id] = atdict

        self.logger.info("collected %d AT reports to post" % len(atreports))

        # If we do not want to submit anything, or if there's nothing to submit
        if len(atreports) == 0 or (not self.submit_tns):
            self.logger.info(
                "submit_tns config parameter is False or there's nothing to submit",
                extra={
                    "n_reports": len(atreports),
                    "submit_tns": self.submit_tns,
                },
            )
            return

        # atreports is now a dict with tran_id as keys and atreport as keys
        # what we need is a list of dicts with form {'at_report':atreport }
        # where an atreport is a dictionary with increasing integer as keys and atreports as values
        atreportlist = [
            {
                "at_report": {
                    i: report
                    for chunk in chunks(atreports.values(), 1)
                    for i, report in enumerate(chunk)
                }
            }
        ]
        tnsreplies = self.client.sendReports(atreportlist)

        # Now go and check and create journal updates for the cases where SN was added
        for tran_id in atreports.keys():
            ztf_name = to_ztf_id(tran_id)
            if not ztf_name in tnsreplies.keys():
                self.logger.info("No TNS add reply", extra={"tranId": tran_id})
                continue


            # Create new journal entry assuming we submitted or found a name
            if "TNSName" in tnsreplies[ztf_name][1].keys():
                gen.send((
                    tran_id,
                    StockAttributes(
                        journal=JournalAttributes(
                            extra={
                                "tnsName": tnsreplies[ztf_name][1]["TNSName"],
                                "tnsInternal": ztf_name,
                                "tnsSubmitresult": tnsreplies[ztf_name][0],
                                "tnsSender": self.tns_api_key.get()["name"],
                            },
                        ),
                        tag="TNS_SUBMITTED",
                        name=tnsreplies[ztf_name][1]["TNSName"],
                    )
                ))


    def report_to_slack(self, atreports: dict[StockId, dict[str, Any]]) -> None:
        self.logger.info("done running T3")

        if not atreports:
            self.logger.info("No atreports collected.")
            return
        elif self.slack_token is None:
            return

        # TODO: to help debugging and verification, we post the collected atreports
        # to the slack, so that we can compare them with what JNo script is doing
        # ALL THE CONTENT OF THIS METHOD SHOULD GO AWAY AS SOON AS WE TRUST THIS T3
        self.logger.warn(
            "Posting collected ATreports to Slack. I'm still running as a test!"
        )

        import datetime
        import io
        import json

        from slack import WebClient
        from slack.errors import SlackClientError
        from slack.web.slack_reponse import SlackResponse

        sc = WebClient(token=self.slack_token.get())

        tstamp = datetime.datetime.today().strftime("%Y-%m-%d-%X")
        atlist = list(atreports.values())
        last = 0
        for ic, atrep in enumerate(chunks(atlist, self.max_slackmsg_size)):

            # add the atreport to a file
            self.logger.info("Posting chunk #%d" % ic)
            filename = "TNSTalker_DEBUG_%s_chunk%d.json" % (tstamp, ic)
            fbuffer = io.StringIO(filename)
            json.dump(atrep, fbuffer, indent=2)

            # upload the file with the at reports
            first = last
            last += len(atrep)
            msg = (
                "A total of %d atreports found by TNSTalker T3. Here's chunk #%d (reports from %d to %d)"
                % (len(atreports), ic, first, last)
            )
            api = sc.files_upload(
                channels=[self.slack_channel],
                title="TNSTalker_%s_chunk%d" % (tstamp, ic),
                initial_comment=msg,
                username=self.slack_username,
                as_user=False,
                filename=filename,
                filetype="javascript",
                file=fbuffer.getvalue(),
            )
            assert isinstance(api, SlackResponse)
            if not api["ok"]:
                raise SlackClientError(api["error"])

        self.logger.warn(
            f"DONE DEBUG Slack posting. Look at {self.slack_channel} for the results"
        )
