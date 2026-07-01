from collections.abc import Mapping, Sequence
from typing import Any

from ampel.contrib.hu.model.LSSTReport import Feature, Host
from ampel.contrib.hu.t2.T2LSSTReport import T2LSSTReport
from ampel.contrib.hu.t2.T2NuclearFilter import T2NuclearFilter


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

        template = Feature(
            name="template_flux",
            version=T2NuclearFilter.version,
            features=nuclear_filter["template_flux"],
        )

        mean_position = Feature(
            name="mean_position",
            version=T2NuclearFilter.version,
            features={k: nuclear_filter[k] for k in ["mean_ra", "mean_dec"]},
        )

        return [host, template, mean_position]
