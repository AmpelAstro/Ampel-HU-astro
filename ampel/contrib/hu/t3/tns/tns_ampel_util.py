#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/tns/tns_ampel_util.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                1.03.2024
# Last Modified Date:  1.03.2024
# Last Modified By:    jnordin@physik.hu-berlin.de

# Methods for converting AMPEL content to TNS compatible data.

from collections.abc import Sequence
from typing import Any

import numpy as np

from ampel.content.DataPoint import DataPoint
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper

TNSFILTERID = {1: "110", 2: "111", 3: "112"}
AT_REPORT_FORM = "bulk-report"
AT_REPORT_REPLY = "bulk-report-reply"
TNS_ARCHIVE = {"OTHER": "0", "SDSS": "1", "DSS": "2"}

# Default atdict settings for ZTF
ZTF_TNS_AT: dict = {  # Default values to tag ZTF detections / ulims
    "flux_units": "1",
    "instrument_value": "196",
    "exptime": "30",
    "Observer": "Robot",
}


def ztfdps_to_tnsdict(
    dps: Sequence[DataPoint] | None,
    max_maglim: float = 19.5,
) -> None | dict[str, Any]:
    """
    Collect ZTF data needed for the atreport. Return None in case
    you have to skip this transient for some reason.
    """

    if not dps:
        return None

    zdps = [dp for dp in dps if "ZTF" in dp["tag"] and "magpsf" in dp["body"]]
    zuls = [
        dp
        for dp in dps
        if "ZTF" in dp["tag"]
        and "magpsf" not in dp["body"]
        and dp["body"]["diffmaglim"] >= max_maglim
    ]

    if len(zdps) == 0:
        return None

    ra = np.mean([dp["body"]["ra"] for dp in zdps])
    dec = np.mean([dp["body"]["dec"] for dp in zdps])

    names: list[str] = []
    for dp in zdps:
        if isinstance(dp["stock"], int):
            names.append(ZTFIdMapper.to_ext_id(dp["stock"]))
        elif isinstance(dp["stock"], Sequence):
            names.extend([ZTFIdMapper.to_ext_id(stock) for stock in dp["stock"]])
    ztfnames = set(names)


    # Start defining AT dict: name and position
    atdict: dict[str, Any] = {}
    atdict["ra"] = {"value": ra, "error": 1.0, "units": "arcsec"}
    atdict["dec"] = {"value": dec, "error": 1.0, "units": "arcsec"}
    atdict["internal_name"] = next(iter(ztfnames))

    # Add information on the latest SIGNIFICANT non detection.
    last_non_obs = 0
    if len(zuls) > 0:
        last_ulim = sorted(zuls, key=lambda x: x["body"]["jd"])[-1]
        last_non_obs = last_ulim["body"]["jd"]
        atdict["non_detection"] = {
            "obsdate": last_ulim["body"]["jd"],
            "limiting_flux": last_ulim["body"]["diffmaglim"],
            "filter_value": TNSFILTERID.get(last_ulim["body"]["fid"]),
        }
    else:
        atdict["non_detection"] = {
            "archiveid": "0",
            "archival_remarks": "ZTF non-detection limits not available",
        }

    atdict["non_detection"].update(ZTF_TNS_AT)  # Add the default ZTF values

    # now add info on photometric detections: consider only candidates which
    # have some consecutive detection after the last ulim
    atdict["photometry"] = {"photometry_group": {}}
    atdict["discovery_datetime"] = 10**30
    for dp in zdps:
        if dp["body"]["jd"] < last_non_obs:
            continue

        # Lets create a few photometry points
        # Note: previously had a cap on the number of dps that could be included. *should* be unnecessary.
        photdict = {
            "obsdate": dp["body"]["jd"],
            "flux": float("{0:.2f}".format(dp["body"]["magpsf"])),  # noqa: UP030
            "flux_error": float("{0:.2f}".format(dp["body"]["sigmapsf"])),  # noqa: UP030
            "limiting_flux": float("{0:.2f}".format(dp["body"]["diffmaglim"])),  # noqa: UP030
            "filter_value": TNSFILTERID.get(dp["body"]["fid"]),
        }
        if dp["body"]["jd"] < atdict["discovery_datetime"]:
            atdict["discovery_datetime"] = dp["body"]["jd"]
        photdict.update(ZTF_TNS_AT)
        atdict["photometry"]["photometry_group"][
            len(atdict["photometry"]["photometry_group"])
        ] = photdict

    return atdict


def get_tns_t2remarks(tview: TransientView) -> None | dict[str, Any]:
    """
    Inspect t2results, and extract TNS remarks when warranted.
    """
    # Tag things close to SDSS nuclei
    nuclear_dist = 1.0

    # Start building dict with remarks
    remarks: dict[str, Any] = {"remarks": ""}

    # Ampel Z
    t2res = tview.get_t2_body(unit="T2DigestRedshifts")
    if isinstance(t2res, dict) and t2res.get("ampel_z", -10) > 0:
        remarks["remarks"] = remarks["remarks"] + "AmpelZ{:.3f} (N{}) ".format(
            t2res["ampel_z"], t2res["group_z_nbr"]
        )

    # T2CatalogMatch
    cat_res = tview.get_t2_body(unit="T2CatalogMatch")
    if isinstance(cat_res, dict):
        # Check redshift
        nedz = cat_res.get("NEDz", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if sdss_spec:
            remarks["remarks"] = (
                remarks["remarks"] + "SDSS spec-z %.3f. " % (sdss_spec["z"])
            )
        elif nedz:
            remarks["remarks"] = remarks["remarks"] + "NED z %.3f. " % (nedz["z"])

        # tag AGNs
        milliquas = cat_res.get("milliquas", False)
        if (
            milliquas
            and milliquas["redshift"] is not None
            and milliquas["redshift"] > 0
        ) or (sdss_spec and sdss_spec["bptclass"] in [4, 5]):
            remarks["remarks"] = (
                remarks["remarks"] + "Known SDSS and/or MILLIQUAS QSO/AGN. "
            )
            remarks["at_type"] = 3

        # tag nuclear
        sdss_dr10 = cat_res.get("SDSSDR10", False)
        if (
            sdss_dr10
            and sdss_dr10["type"] == 3
            and sdss_dr10["dist2transient"] < nuclear_dist
        ):
            remarks["remarks"] = (
                remarks["remarks"] + "Close to core of SDSS DR10 galaxy. "
            )
            remarks["at_type"] = 4

        # Note: removed the tag of noisy gaia data (check T2TNSEval)

    if len(remarks["remarks"]) == 0:
        return None

    return remarks
