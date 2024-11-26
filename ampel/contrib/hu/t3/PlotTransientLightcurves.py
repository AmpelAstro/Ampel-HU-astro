#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/PlotTransientLightcurves.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                11.06.2018
# Last Modified Date:  30.07.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

import gzip
import io
import os
from collections.abc import Generator, Iterable
from contextlib import nullcontext
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy import visualization
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
from slack_sdk.web import WebClient
from ztfquery.utils.stamps import get_ps_stamp

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import T3Send, UBson
from ampel.view.TransientView import TransientView


# Base structure from nuztf (Stein, Reusch)
def fig_from_fluxtable(
    name: str,
    ampelid: str,
    ra: float,
    dec: float,
    fluxtable: Table,
    ztfulims: Table | None = None,
    figsize: tuple = (8, 5),
    title: None | str = None,
    tnsname: str | None = None,
    fritzlink: bool = True,
    attributes: list = [],  # noqa: B006
    cutouts: dict | None = None,
    mag_range: None | list = None,
    z: float | None = None,
    zp: float = 25.0,
    legend: bool = False,
    grid_interval: None | int = None,
    t_0_jd: None | float = None,
    cutout_cache_dir: None | str = ".",
):
    """
    Create a lightcurve figure (in mag space) based on
    flux tables from AMPEL
    """

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Transform to mag (option to plot flux?). Will fail for negative values...
    fluxtable["mag"] = -2.5 * np.log10(fluxtable["flux"]) + zp
    fluxtable["magerr"] = np.abs(
        -2.5 * fluxtable["fluxerr"] / (fluxtable["flux"] * np.log(10))
    )

    if z is not None and np.isnan(z):
        z = None

    # Keys as provided by tabulators
    BANDPASSES = {
        "ztfg": {"label": "ZTF g", "c": "green"},
        "ztfr": {"label": "ZTF R", "c": "red"},
        "ztfi": {"label": "ZTF i", "c": "orange"},
    }

    fig = plt.figure(figsize=figsize)

    # Prepare plot sections
    if cutouts:
        lc_ax1 = fig.add_subplot(5, 4, (9, 19))
        cutoutsci = fig.add_subplot(5, 4, (1, 5))
        cutouttemp = fig.add_subplot(5, 4, (2, 6))
        cutoutdiff = fig.add_subplot(5, 4, (3, 7))
        cutoutps1 = fig.add_subplot(5, 4, (4, 8))
    else:
        lc_ax1 = fig.add_subplot(1, 1, 1)
        fig.subplots_adjust(top=0.8, bottom=0.15)

    plt.subplots_adjust(wspace=0.4, hspace=1.8)

    if cutouts:
        for ax_, type_ in zip(
            [cutoutsci, cutouttemp, cutoutdiff],
            ["Science", "Template", "Difference"],
            strict=False,
        ):
            create_stamp_plot(cutouts=cutouts, ax=ax_, cutout_type=type_)

        if cutout_cache_dir is not None:
            img_cache = os.path.join(cutout_cache_dir, f"{name}_PS1.png")

            if not os.path.isfile(img_cache):
                img = get_ps_stamp(ra, dec, size=240, color=["y", "g", "i"])
                img.save(img_cache)
            else:
                from PIL import Image

                img = Image.open(img_cache)
        else:
            img = get_ps_stamp(ra, dec, size=240, color=["y", "g", "i"])

        cutoutps1.imshow(np.asarray(img))
        cutoutps1.set_title("PS1", fontdict={"fontsize": "small"})
        cutoutps1.set_xticks([])
        cutoutps1.set_yticks([])

    # If redshift is given, calculate absolute magnitude via luminosity distance
    # and plot as right axis
    if z is not None:
        dist_l = cosmo.luminosity_distance(z).to(u.pc).value

        def mag_to_absmag(mag):
            return mag - 5 * (np.log10(dist_l) - 1)

        def absmag_to_mag(absmag):
            return absmag + 5 * (np.log10(dist_l) - 1)

        lc_ax3 = lc_ax1.secondary_yaxis(
            "right", functions=(mag_to_absmag, absmag_to_mag)
        )

        if not cutouts:
            lc_ax3.set_ylabel("Absolute Magnitude [AB]")

    # Give the figure a title
    if not cutouts:
        if title is None:
            fig.suptitle(f"{name}", fontweight="bold")
        else:
            fig.suptitle(title, fontweight="bold")

    if grid_interval is not None:
        lc_ax1.xaxis.set_major_locator(MultipleLocator(grid_interval))

    lc_ax1.grid(visible=True, axis="both", alpha=0.5)
    lc_ax1.set_ylabel("Magnitude [AB]")

    if not cutouts:
        lc_ax1.set_xlabel("JD")

    # Determine magnitude limits
    if mag_range is None:
        max_mag = np.max(fluxtable["mag"]) + 0.3
        min_mag = np.min(fluxtable["mag"]) - 0.3
        lc_ax1.set_ylim((max_mag, min_mag))
    else:
        lc_ax1.set_ylim((np.max(mag_range), np.min(mag_range)))

    for fid in BANDPASSES:
        # Plot older datapoints
        tempTab = fluxtable[fluxtable["band"] == fid]
        lc_ax1.errorbar(
            tempTab["time"],
            tempTab["mag"],
            tempTab["magerr"],
            color=BANDPASSES[fid]["c"],
            fmt=".",
            label=BANDPASSES[fid]["label"],
            mec="black",
            mew=0.5,
        )

        # Plot upper limits
        if ztfulims is not None:
            tempTab = ztfulims[ztfulims["band"] == fid]
            lc_ax1.scatter(
                tempTab["time"],
                tempTab["diffmaglim"],
                c=BANDPASSES[fid]["c"],
                marker="v",
                s=20.0,
                alpha=0.5,
            )

    if legend:
        plt.legend()

    # Now we create an infobox
    if cutouts:
        info = []

        info.append(name)
        info.append(f"RA: {ra:.8f}")
        info.append(f"Dec: {dec:.8f}")
        info.append("------------------------")

        fig.text(0.77, 0.55, "\n".join(info), va="top", fontsize="medium", alpha=0.5)

    # Add annotations
    # Frits
    if fritzlink:
        lc_ax1.annotate(
            "See On Fritz",
            xy=(0.5, 1),
            xytext=(0.78, 0.10),
            xycoords="figure fraction",
            verticalalignment="top",
            color="royalblue",
            url=f"https://fritz.science/source/{name}",
            fontsize=12,
            bbox=dict(boxstyle="round", fc="cornflowerblue", ec="royalblue", alpha=0.4),
        )
    # TNS
    if tnsname is not None:
        lc_ax1.annotate(
            "See On TNS",
            xy=(0.5, 1),
            xytext=(0.78, 0.05),
            xycoords="figure fraction",
            verticalalignment="top",
            color="royalblue",
            url=f"https://www.wis-tns.org/object/{tnsname}",
            fontsize=12,
            bbox=dict(boxstyle="round", fc="cornflowerblue", ec="royalblue", alpha=0.4),
        )

    # Catalog info through T2CatalogMatch results. Add to info above, or keep some of it? At least check TNS somehow
    if len(attributes) > 0:
        ypos = 0.975 if cutouts else 0.035
        fig.text(
            0.5,
            ypos,
            " - ".join(attributes),
            va="top",
            ha="center",
            fontsize="medium",
            alpha=0.5,
        )

    if t_0_jd is not None:
        lc_ax1.axvline(t_0_jd, linestyle=":")
    else:
        t_0_jd = np.mean(fluxtable["time"])

    # Ugly hack because secondary_axis does not work with astropy.time.Time datetime conversion
    jd_min = min(np.min(fluxtable["time"]), t_0_jd)
    if ztfulims is not None:
        jd_min = min(np.min(ztfulims["time"]), jd_min)
    jd_max = max(np.max(fluxtable["time"]), t_0_jd)
    length = jd_max - jd_min

    lc_ax1.set_xlim((jd_min - (length / 20), jd_max + (length / 20)))

    lc_ax2 = lc_ax1.twiny()

    lc_ax2.scatter(  # type: ignore
        [Time(x, format="jd").datetime for x in [jd_min, jd_max]], [20, 20], alpha=0
    )
    lc_ax2.tick_params(axis="both", which="major", labelsize=6, rotation=45)
    lc_ax1.tick_params(axis="x", which="major", labelsize=6, rotation=45)
    lc_ax1.ticklabel_format(axis="x", style="plain")
    lc_ax1.tick_params(axis="y", which="major", labelsize=9)

    if z is not None:
        lc_ax3.tick_params(axis="both", which="major", labelsize=9)

    axes = [lc_ax1, lc_ax2, lc_ax3] if z is not None else [lc_ax1, lc_ax2]

    return fig, axes


def create_stamp_plot(cutouts: dict, ax, cutout_type: str):
    """
    Helper function to create cutout subplot.
    Cutouts assumed to be a dict of the type returned from
    the ZTFCutoutImages complement:
      {'candid': {'cutoutScience': b'..', 'cutoutTemplate': ...} }

    Grabbing images of the first candid available.

    cutout_type assumed to be one of Science, Template, Difference
    """

    data = next(iter(cutouts.values()))[f"cutout{cutout_type}"]

    with gzip.open(io.BytesIO(data), "rb") as f:
        data = fits.open(io.BytesIO(f.read()), ignore_missing_simple=True)[0].data
    vmin, vmax = np.percentile(data[data == data], [0, 100])  # noqa: PLR0124
    data_ = visualization.AsinhStretch()((data - vmin) / (vmax - vmin))
    ax.imshow(
        data_,
        norm=Normalize(*np.percentile(data_[data_ == data_], [0.5, 99.5])),  # noqa: PLR0124
        aspect="auto",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(cutout_type, fontdict={"fontsize": "small"})


# Should be added to ZTFT2Tabulator?
def get_upperlimit_table(
    #        self,
    dps: Iterable[DataPoint] | None,
) -> Table | None:
    if not dps:
        return None

    def filter_limits(dps: Iterable[DataPoint]) -> list[dict]:
        return [
            dp["body"]
            for dp in dps
            if dp["id"] < 0 and "ZTF" in dp["tag"] and "diffmaglim" in dp["body"]
        ]

    ZTF_BANDPASSES = {
        1: {"name": "ztfg"},
        2: {"name": "ztfr"},
        3: {"name": "ztfi"},
    }

    dps_subset = filter_limits(dps)

    if len(dps_subset) == 0:
        return None

    filter_names = [ZTF_BANDPASSES[dp["fid"]]["name"] for dp in dps_subset]
    # fluxlim = np.asarray([10 ** (-((dp["diffmaglim"]) - 25) / 2.5) for dp in dps_subset])

    tab = Table(
        {
            "time": [dp["jd"] for dp in dps_subset],
            "diffmaglim": [dp["diffmaglim"] for dp in dps_subset],
            "band": [ZTF_BANDPASSES[dp["fid"]]["name"] for dp in dps_subset],
            "fluxlim": np.asarray(
                [10 ** (-((dp["diffmaglim"]) - 25) / 2.5) for dp in dps_subset]
            ),
            "zp": [25] * len(filter_names),
            "zpsys": ["ab"] * len(filter_names),
        },
        dtype=("float64", "float64", "str", "float64", "float", "str"),
    )

    return tab  # noqa: RET504


class PlotTransientLightcurves(AbsPhotoT3Unit, AbsTabulatedT2Unit):
    """

    Create a (pdf) plot summarizing lightcurves of candidates provided to the unit.
    Features:
    - Include thumbnails (if provided through the ZTFCutoutImages T3 complement.
    - Include link to TNS (if match existing and provided through TNSNames T3 complement.
    - Upload to slack channel (if token provided)


    """

    # Default path is to create a multi-page pdf
    pdf_path: None | str = None  # Will create random if not set
    titleprefix: str = "AMPEL: "
    # Optionally, save a {stock}.png image of each individual event
    save_png: bool = False
    # Dir for saving png (thumbnails + single event if chosen)
    image_cache_dir: str | None = "./image_cache"

    save_dir: str = "./images"  # dir for saving plots, pdf

    # Should ZTF cutouts be retrieved (requires remote archive access)
    include_cutouts: bool = False

    # Add Fritz link to plot
    fritzlink: bool = True

    # Will post result to Slack channel if a slack channel and a NamedSecret containig the corresponding token is given
    slack_channel: str | None = None
    slack_token: NamedSecret[str] | None = None

    def post_init(self) -> None:
        os.makedirs(self.save_dir, exist_ok=True)
        if self.image_cache_dir:
            os.makedirs(self.image_cache_dir, exist_ok=True)
        # Create temporary path if not set
        if not self.pdf_path:
            import tempfile

            self.pdf_path = tempfile.mkstemp(".pdf", "candidates", self.save_dir)[1]

        # Possibly create a slack client
        if self.slack_channel and self.slack_token is not None:
            self.webclient = WebClient(self.slack_token.get())

    def attributes_from_t2(
        self,
        tview: TransientView,
        nearby_z: float = 0.02,
        snia_minprob: float = 0.7,
        min_kilonovaness=0,
    ) -> tuple[list, Any]:
        """
        Collect information from potential T2 documents,
        return as list of str.
        Partially copied from AstroColibriPublisher. TODO: Join as util function.
        Redshift gets a special treatment, since its ... speical
        """

        attributes = []
        z = None
        # Nearby attribute
        t2res = tview.get_t2_body(unit="T2DigestRedshifts")
        if isinstance(t2res, dict) and t2res.get("ampel_z", -10) > 0:
            attributes.append(
                "AmpelZ{:.3f} N{}".format(t2res["ampel_z"], t2res["group_z_nbr"])
            )
            z = t2res["ampel_z"]
            if t2res.get("ampel_z", 999) < nearby_z:
                attributes.append("Nearby")
        # Infant attribute
        t2res = tview.get_t2_body(unit="T2InfantCatalogEval")
        if isinstance(t2res, dict) and t2res.get("action", False):
            attributes.append("InfantEval")
        # SNIa
        t2res = tview.get_t2_body(unit="T2RunParsnip")
        if (
            isinstance(t2res, dict)
            and "classification" in t2res
            and t2res["classification"]["SNIa"] > snia_minprob
        ):
            attributes.append("ProbSNIa")
        # Kilonovaness
        t2res = tview.get_t2_body(unit="T2KilonovaEval")
        if (
            isinstance(t2res, dict)
            and t2res.get("kilonovaness", -99) > min_kilonovaness
        ):
            attributes.append("Kilonovaness{}".format(t2res["kilonovaness"]))
            attributes.append("LVKmap{}".format(t2res["map_name"]))
        # Linearfit
        t2res = tview.get_t2_body(unit="T2LineFit")
        if isinstance(t2res, dict) and t2res.get("chi2dof", None) is not None:
            attributes.append("LinearChi/dof{:.2}".format(t2res["chi2dof"]))

        t2res = tview.get_t2_body(unit="T2KilonovaStats")
        if isinstance(t2res, dict):
            attributes.append(f"PercentHigher{t2res['gaus_percent']:.5f} ")
            attributes.append(
                f"ExpectedCands{t2res['exp_kn']:.1f}+{t2res['exp_kn_pls']:.1f}-{t2res['exp_kn_min']:.1f} "
            )
            attributes.append(f"DistanceRange{t2res['dist_range']}")
        return (attributes, z)

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None
    ) -> UBson | UnitResult:
        buffer = io.BytesIO()
        with PdfPages(self.pdf_path or buffer) as pdf:
            for tran_view in gen:
                if not tran_view.get_photopoints():
                    self.logger.debug("No photopoints", extra={"stock": tran_view.id})
                    continue
                sncosmo_table = self.get_flux_table(tran_view.get_photopoints())  # type: ignore

                # Collect information
                ampelid = tran_view.id
                (ra, dec) = self.get_pos(tran_view.get_photopoints())  # type: ignore
                name = " ".join(
                    map(str, self.get_stock_name(tran_view.get_photopoints()))  # type: ignore
                )

                # Upper limits (from ZTF)
                # Could immediately subselect to ZTF limits, but keep like this to shift into tabulators
                ulim_table = get_upperlimit_table(tran_view.get_upperlimits())

                # Title
                title = f"{self.titleprefix}: {name!r}-{ampelid!r} @ RA {ra:.3f} Dec {dec:.3f}"

                # Gatter attributes from potential T2 documents
                (attributes, z) = self.attributes_from_t2(tran_view)

                # Check if ZTF name exists in TNS mirror archive
                tnsname = None
                if (
                    isinstance(tran_view.extra, dict)
                    and "TNSReports" in tran_view.extra
                ):
                    tnsname = next(iter(tran_view.extra["TNSReports"])).get(
                        "objname", None
                    )

                if (
                    self.include_cutouts
                    and tran_view.extra
                    and (cutouts := tran_view.extra.get("ZTFCutoutImages", None))
                    is not None
                ):
                    # Complement cutouts worked
                    pass
                else:
                    cutouts = None

                # Create plot
                fig, axes = fig_from_fluxtable(
                    name,
                    str(ampelid),
                    ra,
                    dec,
                    sncosmo_table,
                    ulim_table,
                    title=title,
                    attributes=attributes,
                    fritzlink=self.fritzlink,
                    tnsname=tnsname,
                    z=z,
                    cutouts=cutouts,
                    cutout_cache_dir=self.image_cache_dir,
                )
                pdf.savefig()
                if self.save_png and self.image_cache_dir:
                    plt.savefig(
                        os.path.join(self.image_cache_dir, str(tran_view.id) + ".png")
                    )
                plt.close()

        # Post to slack
        if self.slack_channel is not None and self.slack_token is not None:
            buffer.seek(0)
            with (
                open(self.pdf_path, "rb") if self.pdf_path else nullcontext(buffer)  # type: ignore[attr-defined]
            ) as file:
                self.webclient.files_upload(
                    file=file,
                    filename=self.pdf_path or "candidates.pdf",
                    channels=self.slack_channel,
                    #                thread_ts=self.ts,
                )

        return None
