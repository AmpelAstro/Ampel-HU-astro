from collections.abc import Mapping, Sequence
from typing import Any

from ampel.contrib.hu.model.LSSTReport import (
    Classification,
    Feature,
    Host,
)
from ampel.contrib.hu.t2.T2LSSTReport import T2LSSTReport


class T2InfantReport(T2LSSTReport):
    """
    Report potential infant transients.
    Todo: should we really require RiseDecline? Same information included in photometry extract for table.
    Keep for now, but look into performance.
    Should be possible to run without T2TabulatorRiseDecline - then just remove from report_t2s and set select map to {}.
    In principle, it might be assumed that the T0 filter should be used to only select infant transients.
    """

    report_t2s: Sequence[str] = ["T2DigestRedshifts", "T2TabulatorRiseDecline"]

    # Only transients with features within these limits will be reported
    risedecline_select_map: dict[str, Sequence[float]] = {
        "t_lc": [0, 3],
    }
    redshift_select_map: dict[str, Sequence[float]] = {
        "group_z_nbr": [0, 4],
        "ampel_z": [0, 0.05],
        "ampel_dist": [1.0, 10.0],
    }

    def process_t2s(
        self, report_views: dict[str, Mapping[str, Any]]
    ) -> Sequence[Classification | Host | Feature] | None:
        """
        Process T2 views to extract information to be propagated.
        Will parse results from T2ParsnipRiseDecline and:
        - Check whether minimal quality criteria are met.
        - Collect sparse feature and host information.
        - Look for Parsnip classification results and propagate a summary.
        """

        # To be propagated, collect information.
        unitresults: list[Classification | Host | Feature] = []

        # Check host info
        body = report_views["T2DigestRedshifts"]
        skip = False
        for field, limits in self.redshift_select_map.items():
            if field in body and limits[0] < body[field] < limits[1]:
                continue
            skip = True
        if skip:
            return None

        # Add host information if available
        unitresults.append(
            Host(
                name="AmpelDigestRedshifts",
                redshift=body["ampel_z"],
                redshift_error=body["group_z_precision"],
                distance=body["ampel_dist"],
                info="AMPELz_group" + str(body["group_z_nbr"]),
                source="AmpelZ",
            )
        )

        # Check RiseDecline info if requested
        if "T2TabulatorRiseDecline" in self.report_t2s:
            body = report_views["T2TabulatorRiseDecline"]
            skip = False
            features = {}
            for field, limits in self.risedecline_select_map.items():
                if field in body and limits[0] < body[field] < limits[1]:
                    features[field] = body[field]
                    continue
                skip = True
            if skip:
                return None
            # Base list of features to add
            features.update(
                {
                    field: body[field]
                    for field in ["ndet", "t_lc", "jd_last", "mag_last", "mag_min"]
                    if field in body
                }
            )

        # Add selection features. Todo: add a selection of other interesting ones.
        unitresults.append(
            Feature(name="RiseDeclineFeatures", version="1", features=features)
        )

        return unitresults
