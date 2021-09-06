#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/ampel_tns.py
# License           : BSD-3-Clause
# Author            : Ken Smith
# Date              : May 2016
# Last Modified Date: Feb 2018
# Last Modified By  : Jakob Nordin

# -----------------------------------------------------------------------------
# A python sample code for sending a bulk report to the TNS.
# Original sample code by Ken Smith (May 2016) modified by Jakob Nordin (Feb 2018)
# -----------------------------------------------------------------------------

import json
import re
import time
from typing import Any, Dict, List, Tuple, Union

from requests.models import Response

from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
from ampel.protocol.LoggerProtocol import LoggerProtocol
from requests_toolbelt.sessions import BaseUrlSession

TNSFILTERID = {1: "110", 2: "111", 3: "112"}
AT_REPORT_FORM = "bulk-report"
AT_REPORT_REPLY = "bulk-report-reply"
TNS_ARCHIVE = {"OTHER": "0", "SDSS": "1", "DSS": "2"}
TNS_BASE_URL_SANDBOX = "https://sandbox.wis-tns.org/api/"
TNS_BASE_URL_REAL = "https://www.wis-tns.org/api/"

httpErrors = {
    304: "Error 304: Not Modified: There was no new data to return.",
    400: "Error 400: Bad Request: The request was invalid. An accompanying error message will explain why.",
    403: "Error 403: Forbidden: The request is understood, but it has been refused. An accompanying error message will explain why",
    404: "Error 404: Not Found: The URI requested is invalid or the resource requested, such as a category, does not exists.",
    500: "Error 500: Internal Server Error: Something is broken.",
    503: "Error 503: Service Unavailable.",
    429: "Error 429: Rate Limit Exceeded.",
}

class TNSSession(BaseUrlSession):
    def __init__(self, token: TNSToken, baseURL: str = TNS_BASE_URL_REAL) -> None:
        self.token = token
        super().__init__(baseURL)
        self.headers["User-Agent"] = (
            "tns_marker"
            + json.dumps(
                {"tns_id": self.token.id, "name": self.token.name, "type": "bot"}
            )
        )

    def post(self, method: str, payload: Union[str, Dict[str,Any]], payload_key="data", **kwargs) -> Response:
        for _ in range(10):
            if (response := super().post(
                method,
                files=[
                    ("api_key", (None, self.token.api_key)),
                    (payload_key, (None, json.dumps(payload)))
                ],
                **kwargs
            )).status_code != 429:
                return response
            # back off according to rate-limit headers (see https://www.wis-tns.org/content/tns-newsfeed#comment-wrapper-26286)
            delay = response.headers["x-cone-rate-limit-reset" if response.url.endswith('search') else "x-rate-limit-reset"]
            time.sleep(int(delay))


class TNSClient:
    """Send Bulk TNS Request."""

    def __init__(self, baseURL, logger: LoggerProtocol, token: TNSToken, options={}):
        """
        :param baseURL: Base URL of the TNS API
        :param options: (Default value = {})
        """

        self.logger = logger
        self.session = TNSSession(token, baseURL)

    def jsonResponse(self, r: Response) -> Dict:
        """
        Send JSON response given requests object. Should be a python dict.

        :param r: requests object - the response we got back from the server
        :return d: json response converted to python dict
        """

        d : Dict[str, Any] = {}
        # What response did we get?
        message = None
        status = r.status_code

        if status != 200:
            try:
                message = httpErrors[status]
            except ValueError:
                message = f"Error {status}: Undocumented error"

        if message is not None:
            self.logger.warn("TNS bulk submit: " + message)
            return d

        # Did we get a JSON object?
        try:
            d = r.json()
        except ValueError as e:
            self.logger.error("TNS bulk submit", exc_info=e)
            d = {}
            return d

        # If so, what error messages if any did we get?
        self.logger.info(json.dumps(d, indent=4, sort_keys=True))

        if "id_code" in d.keys() and "id_message" in d.keys() and d["id_code"] != 200:
            self.logger.info(
                "TNS bulk submit: Bad response: code = %d, error = '%s'"
                % (d["id_code"], d["id_message"])
            )
        return d

    def sendBulkReport(self, options) -> dict:
        """
        Send the JSON TNS request
        :param options: the JSON TNS request
        """
        # The requests.post method needs to receive a "files" entry, not "data".  And the "files"
        # entry needs to be a dictionary of tuples.  The first value of the tuple is None.
        self.logger.info("TNS bulk submit: " + "sending request")
        r = self.session.post(AT_REPORT_FORM, options, timeout=300)
        # Construct the JSON response and return it.
        self.logger.info("TNS bulk submit: " + "got response (or timed out)")
        return self.jsonResponse(r)
    
    def addBulkReport(self, report):
        """
        Send the report to the TNS

        :return reportId: TNS report ID

        """
        logger = self.logger
        reply = self.sendBulkReport(report)

        reportId = None

        if reply:
            try:
                reportId = reply["data"]["report_id"]
                logger.info("TNS bulk submit: successful with ID %s" % (reportId))
            except ValueError:
                logger.error("Empty response. Something went wrong. Is the API Key OK?")
            except KeyError:
                logger.error("Cannot find the data key. Something is wrong.")

        return reportId
    
    def sendReports(self, reports: List[Dict]):
        """
        Based on a lists of reportlists, send to TNS.
        Return results for journal entries
        """

        # Submit to TNS
        MAX_LOOP = 25
        SLEEP = 2

        logger = self.logger
        reportresult = {}
        for atreport in reports:

            # Submit a report
            for _ in range(MAX_LOOP):
                reportid = self.addBulkReport(atreport)
                if reportid:
                    logger.info("TNS report ID %s" % (reportid))
                    break
                time.sleep(SLEEP)
            else:
                logger.info("TNS bulk report failed")
                continue

            # Try to read reply
            for _ in range(MAX_LOOP):
                time.sleep(SLEEP)
                response = self.getBulkReportReply(reportid)
                if isinstance(response, list):
                    break
            else:
                logger.info("TNS Report reading failed")
                continue

            # Check whether request was bad. In this case TNS looks to return a list with dicts
            # of failed objects which does not correspond to the order of input atdicts.
            # In any case, nothing in the submit is posted.
            # Hence only checking first element
            bad_request = None
            for key_atprop in ['ra','decl','discovery_datetime']:
                if key_atprop in response[0].keys():
                    try:
                        bad_request = response[0][key_atprop]["value"]["5"]["message"]
                        break
                    except KeyError:
                        pass
            if bad_request is not None:
                logger.info( bad_request )
                continue




            # Parse reply for evaluation
            for k, v in atreport["at_report"].items():
                if "100" in response[k].keys():
                    logger.info(
                        "TNS Inserted with name %s" % (response[k]["100"]["objname"])
                    )
                    reportresult[v["internal_name"]] = [
                        "TNS inserted",
                        {"TNSName": response[k]["100"]["objname"]},
                    ]
                elif "101" in response[k].keys():
                    logger.info(
                        "Already existing with name %s" % (response[k]["101"]["objname"])
                    )
                    reportresult[v["internal_name"]] = [
                        "TNS pre-existing",
                        {"TNSName": response[k]["101"]["objname"]},
                    ]
                    

        return reportresult

    def _bulkReportReply(self, report_id: str) -> Dict[str, Any]:
        """
        Get the report back from the TNS

        :param options: dict containing the report ID
        :return: dict

        """
        self.logger.info("TNS bulk submit: " + "looking for reply report")
        # every TNS endpoint wraps its arguments in `data`, except bulk-report-reply
        r = self.session.post(AT_REPORT_REPLY, report_id, payload_key="report_id", timeout=300)
        self.logger.info("TNS bulk submit: " + "got report (or timed out)")
        return self.jsonResponse(r)

    def getBulkReportReply(self, reportId):
        """
        Get the TNS response for the specified report ID
        :param tnsApiKey: TNS API Key
        :return request: The original request
        :return response: The TNS response
        """

        logger = self.logger
        reply = self._bulkReportReply(reportId)

        response = None
        # reply should be a dict
        if reply and "id_code" in reply.keys() and reply["id_code"] == 404:
            logger.warn(
                f"TNS bulk submit {reportId}: Unknown report. "
                f"Perhaps the report has not yet been processed."
            )

        if reply and "id_code" in reply.keys() and reply["id_code"] == 200:
            try:
                response = reply["data"]["feedback"]["at_report"]
            except KeyError:
                logger.error("TNS bulk submit: cannot find the response feedback payload.")

        # This is a bad request. Still propagate the response for analysis.
        if reply and "id_code" in reply.keys() and reply["id_code"] == 400:
            try:
                response = reply["data"]["feedback"]["at_report"]
            except KeyError:
                logger.error("TNS bulk submit: cannot find the response feedback payload.")


        logger.info(f"TNS bulk submit: got response {response}")

        return response
    
    def getInternalName(self, tns_name: str):
        """
        formerly tnsInternal
        """

        response = self.session.post(
            "get/object",
            {"objname": tns_name, "photometry": 0, "spectra": 0}
        )
        parsed = self.jsonResponse(response)

        if parsed["data"]["reply"]["internal_names"] is None:
            return [], "No internal TNS name"

        return parsed["data"]["reply"]["internal_names"], "Got internal name response"
    
    def search(self, ra: float, dec: float, matchradius: float=5.) -> Tuple[List[str], str]:
        """
        formerly tnsName
        """
        r = self.session.post(
            "get/search",
            {"ra": ra, "dec": dec, "radius": matchradius, "units": "arcsec"}
        )
        parsed = self.jsonResponse(r)

        try:
            tnsnames = [v["prefix"] + v["objname"] for v in parsed["data"]["reply"]]
        except KeyError:
            return [], "Error: No TNS names in response"

        return tnsnames, "Found TNS name(s)"
    
    def getNames(self, ra: float, dec: float, matchradius: float=5.):
        """
        Get names of the first TNS object at location

        formerly get_tnsname
        """
        logger = self.logger
        # Look for TNS name at the coordinate of the transient
        tnsnames, runstatus = self.search(ra, dec, matchradius)
        if re.match("Error", runstatus):
            logger.info("TNS get error", extra={"tns_request": runstatus})
            return None, []
        if len(tnsnames) >= 1:
            tns_name = tnsnames[0]
            if len(tnsnames) > 1:
                logger.debug("Multipe TNS names, choosing first", extra={"tns_names": tnsnames})
        else:
            # No TNS name, then no need to look for internals
            return None, None
        logger.info("TNS get cand id", extra={"tns_name": tns_name})

        # Look for internal name (note that we skip the prefix)
        internal_names, *_ = self.getInternalName(tns_name[2:])
        
        return tns_name, internal_names
