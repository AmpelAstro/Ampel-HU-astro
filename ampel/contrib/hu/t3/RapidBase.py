#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/RapidBase.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 15.07.2019
# Last Modified Date: 06.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import Any, Dict, List, Optional, Tuple, Generator, Union

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.struct.UnitResult import UnitResult
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id
from ampel.types import UBson


# get the science records for the catalog match
def get_catalogmatch_srecs(tran_view, logger):
    cat_res = tran_view.get_science_records(t2_class_name="CATALOGMATCH")
    if len(cat_res) == 0 or cat_res is None or cat_res[-1].get_results() is None:
        logger.info("NO CATALOG MATCH FOR THIS TRANSIENT")
        return {}
    return cat_res[-1].get_results()[-1]["output"]


class RapidBase(AbsPhotoT3Unit):
    """
    Trigger rapid reactions. Intended as base class where the react method can be
    implemented as wished and a testreact method posts test reactions to Slack.

    All transients provided to this unit will trigger reactions. It is assumed that
    selection and filtering has taken place in the T2 and through a
    T3FilteringStockSelector-like selection
    """

    # Unless set, no full reaction will be triggered
    do_react: bool

    # If set, will post trigger to slack
    do_testreact: bool
    # Original
    slack_token: Optional[NamedSecret[str]]
    # Hack
    # slack_token_dict: Optional[ Dict[str,any] ] = {'key':'k','value':'v'}
    # from ampel.secret.DictSecretProvider import NamedSecret
    # slack_token: Secret = NamedSecret(**slack_token_dict)

    slack_channel: str = "#ztf_auto"
    slack_username: str = "AMPEL"

    # List of T2 unit names which should be collected for reaction
    t2info_from: List[str] = []


    def post_init(self) -> None:

        self.name = "RapidBase"
        self.logger.info(f"Initialized T3 RapidBase instance {self.name}")

        # feedback
        for k in self.__annotations__:
            self.logger.info(f"Using {k}={getattr(self, k)}")


    def process(self, gen: Generator[TransientView, JournalAttributes, None]) -> Union[UBson, UnitResult]:
        """
        Loop through transients and check for TNS names and/or candidates to submit
        """

        # We will here loop through transients and react individually
        for tv in gen:

            transientinfo = self.collect_info(tv)
            self.logger.info("reacting", extra={"tranId": tv.id})

            # Ok, so we have a transient to react to
            if self.do_react:
                success, jcontent = self.react(tv, transientinfo)
            # Otherwise, test
            elif self.do_testreact:
                test_success, jcontent = self.test_react(tv, transientinfo)
            
            if jcontent:
                gen.send(JournalAttributes(extra=jcontent))

        return None


    def react(
        self, tran_view: TransientView, info: Optional[Dict[str, Any]]
    ) -> Tuple[bool, Optional[Dict[str, Any]]]:
        """
        Replace with react method adopted to particular facility or output
        """

        raise NotImplementedError("No real reaction implemented in RapidBase")
        return self.test_react(tran_view, info)


    def test_react(
        self, tran_view: TransientView, info: Optional[Dict[str, Any]]
    ) -> Tuple[bool, Optional[Dict[str, Any]]]:
        """ Trigger a test slack report """

        success = False

        if not self.slack_token:
            return False, None

        from slack import WebClient
        from slack.errors import SlackClientError
        from slack.web.slack_response import SlackResponse


        sc = WebClient(self.slack_token.get())
        assert isinstance(tran_view.id, int)
        ztf_name = to_ztf_id(tran_view.id)
        msg = "Ampel RapidReact says: Look up %s. Provided info %s" % (
            ztf_name,
            info,
        )
        api = sc.chat_postMessage(
            channel=self.slack_channel,
            text=msg,
            username=self.slack_username,
            as_user=False,
        )
        assert isinstance(api, SlackResponse)
        if not api["ok"]:
            raise SlackClientError(api["error"])
        else:
            success = True

        description = "Sent SLACK msg"
        self.logger.info(description, extra={"channel": self.slack_channel})

        # Document what we did
        jcontent = {"reaction": description, "success": success}

        return success, jcontent


    def collect_info(self, tran_view: TransientView) -> Optional[Dict[str, Any]]:
        """
        Create an information dict from T2 outputs, which can be used by reactors.
        """

        info: Dict[str, Any] = {}

        for t2unit in self.t2info_from:
            t2_result = tran_view.get_latest_t2_body(unit_id=t2unit)
            if isinstance(t2_result, dict):
               info[t2unit] = t2_result
        return info
