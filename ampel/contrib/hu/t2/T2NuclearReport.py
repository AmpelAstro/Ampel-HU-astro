from collections.abc import Mapping, Sequence
from typing import Any

from ampel.contrib.hu.model.LSSTReport import Host
from ampel.contrib.hu.t2.T2LSSTReport import T2LSSTReport


class T2NuclearReport(T2LSSTReport):
    """
    Propagate information from T2NuclearFilter.
    """

    def process_t2s(
        self, report_views: dict[str, Mapping[str, Any]]
    ) -> Sequence[Host] | None:
        """ """

        nuclear_filter = report_views["T2NuclearFilter"]
        digest_redshifts = report_views["T2DigestRedshifts"]

        # skip
        if not nuclear_filter["passed"]:
            return None

        # host types must be str
        type_info = {
            str(cat_name): {
                str(type_key): str(type_value)
                for type_key, type_value in cat_value.items()
            }
            for cat_name, cat_value in nuclear_filter["host_type"].items()
        }

        host = Host(
            name="T2NuclearFilter",
            redshift=digest_redshifts.get("ampel_z"),
            distance=nuclear_filter["host_dist_arcsec"],
            source=nuclear_filter["host_catalogs"],
            info=type_info,
        )

        return [host]
