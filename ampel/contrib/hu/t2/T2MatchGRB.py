from typing import Any, Iterable, Union
from ampel.types import UBson

import requests
import json

from numpy import mean

from ampel.content.DataPoint import DataPoint
from ampel.struct.UnitResult import UnitResult
from ampel.content.T1Document import T1Document
from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.contrib.hu.util.AmpelHealpix import AmpelHealpix

import astropy.time as atime
import astropy.units as aunits
from astropy.coordinates import SkyCoord


class T2MatchGRB(AbsStateT2Unit, AbsTabulatedT2Unit):
    trigger_jd: float or None = 0
    map_dir: str or None = "./"
    map_name: str or None = None

    ra: float = 0
    dec: float = 0
    radius: float = 100000
    before_time: float = 1  # in days
    after_time: float = 3  # in days

    def post_init(self):
        
        # need to get trigger_time of GW event from healpix map
        if self.map_name is not None or self.trigger_jd is None:
            ah = AmpelHealpix(map_name=self.map_name, map_url="", save_dir=self.map_dir)
            # map_hash = ah.process_map()
            self.trigger_jd = ah.get_triggertime()

        # calculate timeframe to check for GRB events
        self.before_iso = (
            atime.Time(self.trigger_jd, format="jd") - aunits.day * self.before_time
        ).iso
        self.after_iso = (
            atime.Time(self.trigger_jd, format="jd") + aunits.day * self.after_time
        ).iso

        print("T2MATCHGRB:: ", self.before_iso, self.after_iso, self.trigger_jd)

        # get GRB events in timeframe
        self.astrocolibri_allsky()

    def process(
        self,
        compound: T1Document,
        datapoints: Iterable[DataPoint],
    ) -> Union[UBson, UnitResult]:
        
        tmp_skycoord = None

        results = { "temporal_grb": [] }

        # get alert coordinates
        for point in datapoints:
            if point["body"].get("ra"):
                tmp_skycoord = SkyCoord(
                    ra=point["body"]["ra"] * aunits.degree,
                    dec=point["body"]["dec"] * aunits.degree,
                    frame="icrs",
                )
                break
        
        if tmp_skycoord is None:
            self.logger.info("No coordinates to compare.")
            return results

        # for all grb events in timeframe, check if alert coordinates overlap within 1 sigma
        for grb_event in self.event_list:
            #print(grb_event)
            tmp_sep = tmp_skycoord.separation(grb_event["skycoord"])
            #print(tmp_sep)
            #if tmp_sep.value <= grb_event["err"]:
            

            append_match = grb_event.copy()
            append_match["separation"] = tmp_sep.value
            append_match.pop("skycoord", None)
            results["temporal_grb"].append(append_match)
        #print("T2MATCHGRB results:: ", results)
        return results

    def astrocolibri_allsky(self) -> None:
        """
        Request and store all GRB events during -1/+3 timeframe of GW event via astrocolibri.
        """

        grb_filter = {
            "type": "FieldSpecification",
            "field": "type",
            "value": "grb",
            "operation": "==",
            "typeField": "string",
        }

        # Base URL of the API
        url = "https://astro-colibri.science/cone_search"

        # Request parameters (headers, body)
        headers = {"Content-Type": "application/json"}
        body = {
            "filter": grb_filter,
            "time_range": {
                "max": self.after_iso,
                "min": self.before_iso,
            },
            "properties": {
                "type": "cone",
                "position": {"ra": self.ra, "dec": self.dec},
                "radius": self.radius,
            },
        }

        # Perform the POST request
        response = requests.post(
            url, headers=headers, data=json.dumps(body), timeout=20
        )

        # Process the response
        if response.status_code == 200:
            self.logger.debug("Astrocolibri: Response successfully received.")
            events = response.json()["voevents"]
            # print('number of events: ' + str(len(events)))
            # print(events)  # only show transient events. You can also access catalog sources
        else:
            self.logger.info(
                "Astrocolibri: Request did NOT succeed : ", response.status_code
            )
            return

        event_list = []
        for event in events:
            if event["source_name"] == "":
                continue
            tmp_dict = {}
            # print(event.keys())
            tmp_dict["ra"] = event["ra"]
            tmp_dict["dec"] = event["dec"]
            tmp_dict["err"] = event["err"]
            tmp_dict["timestamp"] = event["timestamp"]
            tmp_dict["time_jd"] = atime.Time(
                tmp_dict["timestamp"] * 10**-3 * aunits.s, format="unix"
            ).jd
            tmp_dict["comment"] = event["comment"]
            tmp_dict["skycoord"] = SkyCoord(
                ra=tmp_dict["ra"] * aunits.degree,
                dec=tmp_dict["dec"] * aunits.degree,
                frame="icrs",
            )
            tmp_dict["source_name"] = event["source_name"]
            event_list.append(tmp_dict)
        self.event_list = event_list
        #print("T2MATCHGRB:: ", event_list)
