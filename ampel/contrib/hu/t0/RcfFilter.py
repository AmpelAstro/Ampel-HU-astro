#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : RcfFilter.py
# License           : BSD-3-Clause
# Author            : j. nordin <j.nordin@physik.hu-berlin.de>
# Date              : 06.07.2019
# Last Modified Date: 08.11.2022
# Last Modified By  : j. nordin <j.nordin@physik.hu-berlin.de>

import enum
from typing import Any

from astropy.coordinates import SkyCoord

from ampel.abstract.AbsAlertFilter import AbsAlertFilter
from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol


class RcfFilter(AbsAlertFilter):
    """
    Filter for the ZTF Redshift Completeness Factor program.
    """

    min_ndet: int = 1  #: number of previous detections
    min_dist_to_sso: float  #: distance to nearest solar system object [arcsec]
    min_gal_lat: float  #: minium distance from galactic plane. set to negative to disable cut.
    min_age: float  #: Age as calculated based on previous PP in alert
    max_ipac_age: float  #: Age as calculated based on alert keywords
    max_magpsf: float
    min_rb: float = 0.0  #: real bogus score
    min_drb: float = 0.0  #: deep learning real bogus score

    def post_init(self) -> None:
        # To make this tenable we should create this list dynamically depending
        # on what entries are required by the filter. Now deciding not to
        # include drb in this list, eg.
        self.keys_to_check = (
            "fwhm",
            "elong",
            "magdiff",
            "nbad",
            "distpsnr1",
            "sgscore1",
            "distpsnr2",
            "sgscore2",
            "distpsnr3",
            "sgscore3",
            "isdiffpos",
            "ra",
            "dec",
            "rb",
            "drb",
            "ssdistnr",
            "jdstarthist",
        )

    class RejectionCode(enum.IntEnum):
        min_ndet = -1
        has_keys = -2
        isdiffpos = -3
        max_magpsf = -4
        min_gal_lat = -5
        min_ssdistnr = -6
        min_age = -7
        max_ipac_age = -8
        star_under = -9
        is_not_real = -10
        is_bright_star = -11
        is_variable_star = -12

    def _alert_has_keys(self, alert: dict[str, Any]) -> bool:
        """
        check that given photopoint contains all the keys needed to filter
        """
        for el in self.keys_to_check:
            if el not in alert:
                self.logger.debug(None, extra={"missing": el})
                return False
            if alert[el] is None:
                self.logger.debug(None, extra={"isNone": el})
                return False
        return True

    def get_galactic_latitude(self, alert: dict[str, Any]) -> float:
        """
        compute galactic latitude of the transient
        """
        coordinates = SkyCoord(alert["ra"], alert["dec"], unit="deg")
        b = coordinates.galactic.b.deg
        return b

    def previous_pointsource(self, alert: dict[str, Any], age: float) -> bool:
        """
        Sets of cases under which a previous source likely existed.

        Reject if a previous point source (star) exists beneath. Assumed to be the case if any of:
        - sgscore1>0.76 & 0<distpsnr1<2
        - sgscore>0.2 & distpsnr1 < 1 & srmag1 > 0 &
          (szmag1 > 0 & srmag1 - szmag1 > 3.0) or (simag1 > 0 & srmag1 - simag1 > 3.0)
        """

        if (
            alert["sgscore1"] > 0.76
            and 0 < alert["distpsnr1"]
            and alert["distpsnr1"] < 2
        ):
            return True
        if (
            alert["sgscore1"] > 0.2
            and 0 < alert["distpsnr1"]
            and alert["distpsnr1"] < 1
            and alert["srmag1"] > 0
        ):
            if alert["szmag1"] > 0 and (alert["srmag1"] - alert["szmag1"]) > 3.0:
                return True
            if alert["simag1"] > 0 and (alert["srmag1"] - alert["simag1"]) > 3.0:
                return True

        # No clause fulfilled
        return False

    def is_not_real(self, alert: dict[str, Any], age: float) -> bool:
        """
        # Require a high RB score.
        # Generally requires >0.2, but is more stringent if close to a bright catalogued star.
        real = False;
        if (rbscore > 0.2) {
                real = True;
        }
        if (rbscore < 0.35) {
                if (neargaia < 1.0 and neargaia > 0. and maggaia < 17.0 and maggaia > 0.) {real = False;}
                if (distpsnr1 < 1.0 and (srmag1 > 0 and srmag1 < 17.0 or simag1 > 0 and simag1 < 17.0 or szmag1 > 0 and szmag1 < 16.5) and sgscore1 > 0.49) {real = False;}
        }
        if (rbscore < 0.45) {
                if (neargaia < 1.5 and neargaia > 0. and maggaia < 15.5 and maggaia > 0.) {real = False;}
                if (distpsnr1 < 1.5 and (srmag1 > 0 and srmag1 < 15.5 or simag1 > 0 and simag1 < 15.5 or szmag1 > 0 and szmag1 < 15.0) and sgscore1 > 0.49) {real = False;}
        }
        if (drb < 0.5) {
                if (distpsnr1 < 3.0 and ps1mag < 16 and age > 90) { real = False; }
                if (distpsnr1 < 1.1 and ps1mag < 18 and age > 90) { real = False; }
        }
        if (drb < 0.8) {
                if (distpsnr1 < 1.5 and ps1mag < 15.5 and age > 90) { real = False; }
                if (distpsnr1 < 0.8 and ps1mag < 17.5 and age > 90) { real = False; }
        }
        """

        ps1mag = alert["srmag1"]
        if ps1mag <= 0 or ps1mag > 30:
            ps1mag = alert["simag1"]
        if ps1mag <= 0 or ps1mag > 30:
            ps1mag = alert["sgmag1"]
        if ps1mag <= 0 or ps1mag > 30:
            ps1mag = alert["szmag1"]
        if ps1mag <= 0:
            ps1mag = 99

        if alert["rb"] < 0.2:
            return True
        if alert["rb"] < 0.35:
            if (
                0 < alert["neargaia"]
                and alert["neargaia"] < 1.0
                and alert["maggaia"] > 0
                and alert["maggaia"] < 17.0
            ):
                return True
            if (
                alert["distpsnr1"] > 0
                and alert["distpsnr1"] < 1.0
                and alert["sgscore1"] > 0.49
            ):
                if alert["srmag1"] > 0 and alert["srmag1"] < 17:
                    return True
                if alert["simag1"] > 0 and alert["simag1"] < 17:
                    return True
                if alert["szmag1"] > 0 and alert["szmag1"] < 16.5:
                    return True
        if alert["rb"] < 0.45:
            if (
                0 < alert["neargaia"]
                and alert["neargaia"] < 1.5
                and alert["maggaia"] > 0
                and alert["maggaia"] < 15.5
            ):
                return True
            if (
                alert["distpsnr1"] > 0
                and alert["distpsnr1"] < 1.5
                and alert["sgscore1"] > 0.49
            ):
                if alert["srmag1"] > 0 and alert["srmag1"] < 15.5:
                    return True
                if alert["simag1"] > 0 and alert["simag1"] < 15.5:
                    return True
                if alert["szmag1"] > 0 and alert["szmag1"] < 15.0:
                    return True

        if alert["drb"] < 0.5:
            if (
                0 < alert["distpsnr1"]
                and alert["distpsnr1"] < 3.0
                and ps1mag < 16
                and age > 90
            ):
                return True
            if (
                0 < alert["distpsnr1"]
                and alert["distpsnr1"] < 1.1
                and ps1mag < 18
                and age > 90
            ):
                return True

        if alert["drb"] < 0.8:
            if (
                0 < alert["distpsnr1"]
                and alert["distpsnr1"] < 1.5
                and ps1mag < 15.5
                and age > 90
            ):
                return True
            if (
                0 < alert["distpsnr1"]
                and alert["distpsnr1"] < 0.8
                and ps1mag < 17.5
                and age > 90
            ):
                return True

        # Passed tests
        return False

    def is_bright_star(self, alert: dict[str, Any], age: float) -> bool:
        """
        # Could this be residuals from a bright star?
        * Rejected as a bright star if any of:
                - 0 < neargaiabright < 20 & 0 < maggaiabright < 12
                - for i in 1,2,3 and f in r,i,z any of:
                        distpsnri<20 & 0<sfmagi<14 & sgscorei > 0.49
          (original filter only cuts around 10" from z-band and ignored 2,3 for i/z)
        """

        if (
            0 < alert["neargaiabright"]
            and alert["neargaiabright"] < 20
            and 0 < alert["maggaiabright"]
            and alert["maggaiabright"] < 12
        ):
            return True

        for d in ["1", "2", "3"]:
            for f in ["r", "i", "z"]:
                if (
                    0 < alert[f"distpsnr{d}"]
                    and alert[f"distpsnr{d}"] < 20
                    and 0 < alert[f"s{f}mag{d}"]
                    and alert[f"s{f}mag{d}"] < 14
                    and alert[f"sgscore{d}"] > 0.9
                ):
                    return True

        # Passed tests
        return False

    def is_variable_star(
        self, alert: dict[str, Any], age: float, m_peak: float, bright_detections: int
    ) -> bool:
        """
        Rejected as a variable source if any of:
          - age>90 & not peaking below 18.5 & >= 2 prev. det (>=1 if older than 1 yr) & one of:
            - 0<magnr<19.5 & distnr<0.4
            - 0<magnr<17.5 & distnr<0.8
            - 0<magnr<15.5 & distnr<1.2
          - 0<maggaia<17 & 0<neargaia<0.35 & age>30
          - 0<maggaia<19 & 0<neargaia<0.35 & age>300 & magpsf > 18.5
          - 0<maggaia<18 & 0<neargaia<0.2 & age>90
          - 0<magnr<(magpsf-1) & age>90 & distnr<0.5 & not at max
        """

        ps1maxmag = alert["srmag1"]
        if alert["simag1"] > 0 and alert["simag1"] < ps1maxmag:
            ps1maxmag = alert["simag1"]
        if alert["sgmag1"] > 0 and alert["sgmag1"] < ps1maxmag:
            ps1maxmag = alert["sgmag1"]
        if alert["szmag1"] > 0 and alert["szmag1"] < ps1maxmag:
            ps1maxmag = alert["szmag1"]
        if ps1maxmag <= 0:
            ps1maxmag = 99

        if m_peak == alert["magpsf"]:
            is_at_peak = True
        else:
            is_at_peak = False

        if (
            (age > 90 and bright_detections > 2)
            or (age > 365 and bright_detections > 1)
        ) and not (is_at_peak and alert["magpsf"] < 18.5):
            if (
                0 < alert["magnr"]
                and alert["magnr"] < 19.5
                and 0 < alert["distnr"]
                and alert["distnr"] < 0.4
            ):
                return True
            if (
                0 < alert["magnr"]
                and alert["magnr"] < 17.5
                and 0 < alert["distnr"]
                and alert["distnr"] < 0.8
            ):
                return True
            if (
                0 < alert["magnr"]
                and alert["magnr"] < 15.5
                and 0 < alert["distnr"]
                and alert["distnr"] < 1.2
            ):
                return True

        if 0 < alert["maggaia"] and 0 < alert["neargaia"]:
            if alert["neargaia"] < 0.35:
                if age > 30 and alert["maggaia"] < 17:
                    return True
                if age > 300 and alert["maggaia"] < 19 and alert["magpsf"] > 18.5:
                    return True
            if alert["neargaia"] < 0.20:
                if age > 90 and alert["maggaia"] < 18:
                    return True

        if (
            alert["sgscore1"] > 0.25
            and alert["distpsnr1"] < 3
            and age > 90
            and ps1maxmag < 16
        ):
            return True

        if (
            alert["sgscore1"] == 0.5
            and alert["distpsnr1"] < 0.5
            and age > 90
            and ps1maxmag < 17
        ):
            return True

        if (
            age > 90
            and 0 < alert["distnr"]  # shouldn't this be distnr?  was distnbr
            and alert["distnr"] < 0.5
            and not is_at_peak
        ):
            if 0 < alert["magnr"] and alert["magnr"] < (alert["magpsf"] - 1):
                return True

        # Passed tests
        return False

    def process(self, alert: AmpelAlertProtocol) -> None | bool | int:
        # BASE REQUIREMENTS

        codes = self.RejectionCode

        # Number of detections (not in default)
        npp = len(alert.datapoints)
        if npp < self.min_ndet:
            self.logger.debug(None, extra={"nDet": npp})
            return codes.min_ndet

        latest = alert.datapoints[0]
        if not self._alert_has_keys(latest):
            return codes.has_keys

        # require positive detection
        if latest["isdiffpos"] == "f" or latest["isdiffpos"] == "0":
            self.logger.debug(None, extra={"isdiffpos": latest["isdiffpos"]})
            return codes.isdiffpos

        # check brightness
        if latest["magpsf"] > self.max_magpsf:
            self.logger.debug(None, extra={"magpsf": latest["magpsf"]})
            return codes.max_magpsf

        # cut on galactic latitude
        b = self.get_galactic_latitude(latest)
        if abs(b) < self.min_gal_lat:
            self.logger.debug(None, extra={"gal_plane": abs(b)})
            return codes.min_gal_lat

        # check for closeby ss objects
        if 0 <= latest["ssdistnr"] < self.min_dist_to_sso:
            self.logger.debug(None, extra={"ssdistnr": latest["ssdistnr"]})
            return codes.min_ssdistnr

        # find first detection date (with pos. detection)
        jd_first_pps = 10**30
        m_peak = 100
        bright_detections = 0
        for al in alert.datapoints:
            if (
                "isdiffpos" not in al.keys()
                or al["isdiffpos"] == "f"
                or al["isdiffpos"] == "0"
            ):
                continue
            if al["jd"] < jd_first_pps:
                jd_first_pps = al["jd"]
            if al["magpsf"] < m_peak:
                m_peak = al["magpsf"]
            if al["magpsf"] < self.max_magpsf:
                bright_detections += 1

        # print(latest['jd'])
        age = latest["jd"] - latest["jdstarthist"]

        self.logger.debug(
            f"{alert.id} age {age} m now {latest['magpsf']:.2f} m peak {m_peak:.2f} nbr previous bright {bright_detections}",
            extra={"age": age},
        )
        if age < self.min_age:
            self.logger.debug(None, extra={"age": age})
            return codes.min_age

        if self.max_ipac_age > 0:
            try:
                ipac_age = latest["jdendhist"] - latest["jdstarthist"]
                if self.max_ipac_age < ipac_age:
                    self.logger.debug(None, extra={"age": ipac_age})
                    return codes.max_ipac_age
            except KeyError:
                self.logger.debug(
                    "%s No jd end or start alert keywords. Letting through."
                    % (alert.id)
                )

        # SEARCH POINT SOURCE UNDERNEATH
        if self.previous_pointsource(latest, age):
            self.logger.debug("Star under")
            return codes.star_under

        # CHECK IF REAL
        if self.is_not_real(latest, age):
            self.logger.debug("Not real")
            return codes.is_not_real

        # CHECK IF BRIGHT STAR
        if self.is_bright_star(latest, age):
            self.logger.debug("Bright star")
            return codes.is_bright_star

        # CHECK IF VARIABLE STAR
        if self.is_variable_star(latest, age, m_peak, bright_detections):
            self.logger.debug("Variable star")
            return codes.is_variable_star

        # congratulation alert! you made it!
        self.logger.debug("Alert accepted", extra={"latestPpId": latest["candid"]})

        return True
