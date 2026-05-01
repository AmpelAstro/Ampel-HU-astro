#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/PlotNuclearFilterLightcurves.py
# License:             BSD-3-Clause
# Author:              Jannis Necker <jannis.necker@gmail.com>
# Date:                30.04.2026
# Last Modified Date:  30.04.2026
# Last Modified By:    Jannis Necker <jannis.necker@gmail.com>

import base64
import gzip
import io
import os
import tempfile
from collections.abc import Generator, Iterable, Mapping, Sequence
from contextlib import suppress
from gzip import BadGzipFile
from typing import Any

import backoff
import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import units as u
from astropy import visualization
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from requests_toolbelt.sessions import BaseUrlSession
from scipy import ndimage

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.util.catalog_column_info import (
    get_catalog_position_unit_map,
)
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import T3Send, UBson
from ampel.view.TransientView import TransientView


def fig_from_fluxtable(
    name: str,
    ampelid: str,
    ra: float,
    dec: float,
    fluxtable: Table,
    ztfulims: Table | None = None,
    tnsname: str | None = None,
    fritzlink: bool = True,
    attributes: list[str] | None = None,
    nuclear_filter_res: dict[str, Any] | None = None,
    photz_list: list[str] | None = None,
    cutouts: Mapping[str, Mapping[str, bytes]] | None = None,
    mag_range: None | list = None,
    z: float | None = None,
    zp: float = 25.0,
    legend: bool = True,
    title: None | str = None,
    grid_interval: None | int = None,
    t_0_jd: None | float = None,
    finder_cache_dir: None | str = ".",
    cutout_cache_dir: None | str = ".",
    cutout_cache_key: str | None = None,
    stacking: bool = False,
    stacking_window_hours: float = 2.0,
    stacking_alpha: float = 0.22,
    finder_matches: list[dict[str, Any]] | None = None,
):
    """
    Create a lightcurve figure (in mag space) based on flux tables from AMPEL.

    If cutouts are provided, a multi-panel layout is used:
      - Science/Template/Difference cutouts
      - Finder stamp with marker of potential host matches
      - Lightcurve with optional model curve and rolling stacking of data points
      - Optional classification radar + text table
      - General info box (coords, attributes, photo-z)

    If no cutouts are provided, only the lightcurve is shown.
    """
    attributes = attributes or []
    photz_list = photz_list or []

    has_cutouts = bool(cutouts)
    has_matches = bool(finder_matches)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Optional stacking (rolling binning per band)
    raw_fluxtable = fluxtable  # keep for background plotting
    if stacking:
        fluxtable = stack_fluxtable_rolling(
            fluxtable,
            window_hours=stacking_window_hours,
        )

    # Prefer per-point zp from tabulator output, fallback to function argument
    if "zp" in fluxtable.colnames:
        zps = np.asarray(fluxtable["zp"], dtype=float)
    else:
        zps = np.full(len(fluxtable), float(zp), dtype=float)

    # -- Transform flux to mags (AB) --
    flux = np.asarray(fluxtable["flux"], dtype=float)
    fluxerr = np.asarray(fluxtable["fluxerr"], dtype=float)

    fluxtable["mag"] = -2.5 * np.log10(flux) + zps

    # Asymmetric magnitude errors from flux±fluxerr
    flux_plus = flux + fluxerr
    flux_minus = flux - fluxerr

    # Guard against non-positive fluxes (log undefined)
    mag_plus = np.full_like(flux, np.nan, dtype=float)
    mag_minus = np.full_like(flux, np.nan, dtype=float)

    ok_plus = flux_plus > 0
    ok_minus = flux_minus > 0

    mag_plus[ok_plus] = -2.5 * np.log10(flux_plus[ok_plus]) + zps[ok_plus]
    mag_minus[ok_minus] = -2.5 * np.log10(flux_minus[ok_minus]) + zps[ok_minus]

    mag = np.asarray(fluxtable["mag"], dtype=float)

    fluxtable["magerrmin"] = mag - mag_plus
    fluxtable["magerrmax"] = mag_minus - mag

    if z is not None and np.isnan(z):
        z = None

    # Keys as provided by tabulators
    BANDPASSES = {
        # ZTF
        "ztfg": {"label": "ZTF g", "c": "green"},
        "ztfr": {"label": "ZTF r", "c": "red"},
        "ztfi": {"label": "ZTF i", "c": "orange"},
        # LSST
        "lsstu": {"label": "LSST u", "c": "purple"},
        "lsstg": {"label": "LSST g", "c": "blue"},
        "lsstr": {"label": "LSST r", "c": "green"},
        "lssti": {"label": "LSST i", "c": "orange"},
        "lsstz": {"label": "LSST z", "c": "red"},
        "lssty": {"label": "LSST y", "c": "darkred"},
    }

    # General info box
    info: list[str] = [f"ID: {name}", f"RA: {ra:.4f}", f"Dec: {dec:.4f}"]

    if str(ampelid) != str(name):
        info.append("AmpelID:")
        info.append(f"{ampelid}")

    sep = "------------------------"
    info.append(sep)
    if nuclear_filter_res:
        passed = "passed" if nuclear_filter_res["passed"] else "failed"
        info.append(f"Nuclear filter: {passed}")
        if (host_dist := nuclear_filter_res["host_dist_arcsec"]) is not None:
            info.append(f"Host dist: {host_dist:.2f}")
            info.append("Types:")
            for cat_name, cat_info in nuclear_filter_res["host_type"].items():
                cat_info_details = "; ".join(
                    [f"{v} ({k})" for k, v in cat_info.items()]
                )
                if len(cat_info_details):
                    info.append(f"   {cat_name}: {cat_info_details}")
    if attributes:
        info.append(sep)
        info.extend(attributes)
    if photz_list:
        info.append(sep)
        info.append("Photo-z:")
        info.extend(photz_list)

    fig = plt.figure(figsize=(10, 5))
    lc_ax3 = None

    # Helper functions for layout
    def _plot_cutouts_and_get_fov(
        *,
        cutouts,
        axes,
        cache_dir: str | None,
        cache_key: str,
    ) -> float | None:
        """
        Render Science/Template/Difference cutouts into the provided axes and
        return an FOV in arcseconds to use as the finder FOV.
        """
        cutout_fov = None
        for ax_, typ in zip(axes, ["Science", "Template", "Difference"], strict=False):
            fov = create_stamp_plot(
                cutouts=cutouts,
                ax=ax_,
                cutout_type=typ,
                cache_dir=cache_dir,
                cache_key=cache_key,
            )
            if cutout_fov is None and fov is not None and np.isfinite(fov):
                cutout_fov = float(fov)
        return cutout_fov

    def _finder_fov_arcsec(cutout_fov: float | None, factor: float) -> float:
        """
        Compute the finder stamp field-of-view in arcseconds from a cutout FOV.
        Uses a multiplicative factor and clips to a conservative range to avoid extreme values.
        """
        if cutout_fov is None:
            return 45.0
        return float(np.clip(factor * cutout_fov, 10, 1200))

    def _crosshair_gap_frac(cutout_fov: float | None, finder_fov: float) -> float:
        """
        Fractional width of the central crosshair gap so that it matches the
        cutout FoV on the sky.
        """
        if cutout_fov is None or not np.isfinite(cutout_fov) or finder_fov <= 0:
            return 0.2
        return float(np.clip(cutout_fov / finder_fov, 0.05, 0.9))

    def _setup_gs_2x32(fig):
        """
        Create and space the common 2x32 GridSpec used for layouts with cutous.
        """
        gs = GridSpec(
            nrows=2,
            ncols=32,
            figure=fig,
            height_ratios=[1.0, 1.35],
            width_ratios=[1] * 34,
        )
        fig.subplots_adjust(
            left=0.09,
            right=0.98,
            top=0.93,
            bottom=0.10,
            wspace=0.35,
            hspace=0.45,
        )
        return gs

    # merge nuclear filter res into finder mathes
    if nuclear_filter_res and nuclear_filter_res["host_ra"] is not None:
        finder_matches.append(
            {
                "label": "T2NuclearFilter",
                "ra": nuclear_filter_res["host_ra"],
                "dec": nuclear_filter_res["host_dec"],
            }
        )

    # Layout A: Cutouts + finder stamp.
    if has_cutouts:
        assert cutouts is not None
        gs = _setup_gs_2x32(fig)

        # Top row: cutouts + finder + finder legend
        cutoutsci = fig.add_subplot(gs[0, 0:6])
        cutouttemp = fig.add_subplot(gs[0, 6:12])
        cutoutdiff = fig.add_subplot(gs[0, 12:18])
        cutoutfinder = fig.add_subplot(gs[0, 18:24])
        finderleg_ax = fig.add_subplot(gs[0, 24:34]) if has_matches else None

        # Bottom row: lightcurve + attributes text box
        lc_ax1 = fig.add_subplot(gs[1, 0:20])
        attr_ax = fig.add_subplot(gs[1, 24:32])
        attr_ax.axis("off")

        cutout_fov = _plot_cutouts_and_get_fov(
            cutouts=cutouts,
            axes=[cutoutsci, cutouttemp, cutoutdiff],
            cache_dir=cutout_cache_dir,
            cache_key=cutout_cache_key or name,
        )

        finder_fov_arcsec = _finder_fov_arcsec(cutout_fov, factor=1.0)
        crosshair_gap_frac = _crosshair_gap_frac(cutout_fov, finder_fov_arcsec)

        render_finder_stamp(
            cutoutfinder,
            ra,
            dec,
            cache_dir=finder_cache_dir,
            name=name,
            size=240,
            fov_arcsec=finder_fov_arcsec,
            crosshair_gap_frac=crosshair_gap_frac,
            matches=finder_matches,
            legend_ax=finderleg_ax,
        )

        attr_ax.text(
            0.0,
            1.25,
            "\n".join(info),
            va="top",
            ha="left",
            fontsize=8.5,
            alpha=0.6,
            wrap=True,
        )

    # Layout B: Lightcurve + finder + text box (no cutouts).
    else:
        gs = GridSpec(
            nrows=2,
            ncols=32,
            figure=fig,
            height_ratios=[1.0, 1.35],
            width_ratios=[1] * 32,
        )
        fig.subplots_adjust(
            left=0.09,
            right=0.98,
            top=0.83,
            bottom=0.10,
            wspace=0.35,
            hspace=0.45,
        )

        # Left: lightcurve over full height
        lc_ax1 = fig.add_subplot(gs[:, 0:19])

        # Right top: finder + optional legend
        cutoutfinder = fig.add_subplot(gs[0, 23:29])
        finderleg_ax = fig.add_subplot(gs[0, 29:32]) if has_matches else None

        # Right bottom: info / class probabilities
        attr_ax = fig.add_subplot(gs[1, 23:32])
        attr_ax.axis("off")

        render_finder_stamp(
            cutoutfinder,
            ra,
            dec,
            cache_dir=finder_cache_dir,
            name=name,
            size=240,
            fov_arcsec=10,
            crosshair_gap_frac=0.1,
            matches=finder_matches,
            legend_ax=finderleg_ax,
        )

        attr_lines = list(info)

        attr_ax.text(
            0.0,
            1.25,
            "\n".join(attr_lines),
            va="top",
            ha="left",
            fontsize=8.5,
            alpha=0.6,
            wrap=True,
            transform=attr_ax.transAxes,
        )

    # If redshift is given, calculate absolute magnitude via luminosity distance
    if z is not None:
        dist_l = cosmo.luminosity_distance(z).to(u.pc).value

        def mag_to_absmag(mag):
            return mag - 5 * (np.log10(dist_l) - 1)

        def absmag_to_mag(absmag):
            return absmag + 5 * (np.log10(dist_l) - 1)

        lc_ax3 = lc_ax1.secondary_yaxis(
            "right", functions=(mag_to_absmag, absmag_to_mag)
        )
        lc_ax3.set_ylabel("Absolute Magnitude [AB]")

    # Attribute line (kept close to original behavior)
    if attributes:
        fig.text(
            0.5,
            0.985,
            " · ".join(attributes),
            ha="center",
            va="top",
            fontsize=9,
            alpha=0.65,
        )

    if grid_interval is not None:
        lc_ax1.xaxis.set_major_locator(MultipleLocator(grid_interval))

    lc_ax1.grid(visible=True, axis="both", alpha=0.5)
    lc_ax1.set_ylabel("Magnitude [AB]")
    lc_ax1.set_xlabel("JD")

    # Determine magnitude limits
    if mag_range is None:
        max_mag = np.max(fluxtable["mag"]) + 0.3
        min_mag = np.min(fluxtable["mag"]) - 0.3
        lc_ax1.set_ylim((max_mag, min_mag))
    else:
        lc_ax1.set_ylim((np.max(mag_range), np.min(mag_range)))

    # Plot points and upper limits per band
    for fid in BANDPASSES:
        # Background: raw points (if stacking)
        if stacking and raw_fluxtable is not None:
            rawTab = raw_fluxtable[raw_fluxtable["band"] == fid]
            if len(rawTab) > 0:
                # ensure mag columns exist for raw as well
                raw_flux = np.asarray(rawTab["flux"], dtype=float)
                raw_fluxerr = np.asarray(rawTab["fluxerr"], dtype=float)

                if "zp" in rawTab.colnames:
                    raw_zps = np.asarray(rawTab["zp"], dtype=float)
                else:
                    raw_zps = np.full(len(rawTab), float(zp), dtype=float)

                raw_mag = -2.5 * np.log10(raw_flux) + raw_zps

                raw_flux_plus = raw_flux + raw_fluxerr
                raw_flux_minus = raw_flux - raw_fluxerr

                raw_mag_plus = np.full_like(raw_flux, np.nan, dtype=float)
                raw_mag_minus = np.full_like(raw_flux, np.nan, dtype=float)

                okp = raw_flux_plus > 0
                okm = raw_flux_minus > 0

                raw_mag_plus[okp] = -2.5 * np.log10(raw_flux_plus[okp]) + raw_zps[okp]
                raw_mag_minus[okm] = -2.5 * np.log10(raw_flux_minus[okm]) + raw_zps[okm]

                raw_magerrmin = raw_mag - raw_mag_plus
                raw_magerrmax = raw_mag_minus - raw_mag

                lc_ax1.errorbar(
                    rawTab["time"],
                    raw_mag,
                    yerr=[raw_magerrmin, raw_magerrmax],
                    color=BANDPASSES[fid]["c"],
                    fmt=".",
                    markersize=6,
                    alpha=float(stacking_alpha),
                    mec="black",
                    mew=0.3,
                    label=None,  # avoid duplicate legend entries
                    zorder=1,
                )

        # Foreground: (possibly stacked) points
        tempTab = fluxtable[fluxtable["band"] == fid]
        if len(tempTab) > 0:
            lc_ax1.errorbar(
                tempTab["time"],
                tempTab["mag"],
                yerr=[tempTab["magerrmin"], tempTab["magerrmax"]],
                color=BANDPASSES[fid]["c"],
                fmt=".",
                markersize=9,
                label=BANDPASSES[fid]["label"],
                mec="black",
                mew=0.5,
                zorder=2,
            )

        # Upper limits remain unchanged (typically not stacked)
        if ztfulims is not None:
            tempTab = ztfulims[ztfulims["band"] == fid]
            if len(tempTab) > 0:
                lc_ax1.scatter(
                    tempTab["time"],
                    tempTab["diffmaglim"],
                    c=BANDPASSES[fid]["c"],
                    marker="v",
                    s=20.0,
                    alpha=0.5,
                    zorder=0,
                )

    if legend:
        handles, labels = lc_ax1.get_legend_handles_labels()
        if handles:
            lc_ax1.legend(handles, labels, loc="best", fontsize="small")

    # Add annotations
    if fritzlink and not any("lsst" in str(b).lower() for b in fluxtable["band"]):
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

    if t_0_jd is not None:
        lc_ax1.axvline(t_0_jd, linestyle=":")
    else:
        t_0_jd = float(np.mean(fluxtable["time"]))

    # Ugly hack because secondary_axis does not work with astropy.time.Time datetime conversion
    jd_min = min(float(np.min(fluxtable["time"])), float(t_0_jd))
    if ztfulims is not None:
        jd_min = min(float(np.min(ztfulims["time"])), jd_min)
    jd_max = max(float(np.max(fluxtable["time"])), float(t_0_jd))
    length = jd_max - jd_min

    if length > 0:
        lc_ax1.set_xlim((jd_min - (length / 20), jd_max + (length / 20)))
    else:
        lc_ax1.set_xlim((jd_min - 3, jd_max + 3))

    lc_ax2 = lc_ax1.twiny()
    lc_ax2.scatter(
        [Time(x, format="jd").datetime for x in [jd_min, jd_max]], [20, 20], alpha=0
    )

    lc_ax2.tick_params(axis="both", which="major", labelsize=8, rotation=30)
    lc_ax1.tick_params(axis="x", which="major", labelsize=8, rotation=30)
    lc_ax1.ticklabel_format(axis="x", style="plain")
    lc_ax1.tick_params(axis="y", which="major", labelsize=9)

    if lc_ax3 is not None:
        lc_ax3.tick_params(axis="both", which="major", labelsize=9)

    axes = [lc_ax1, lc_ax2, lc_ax3] if lc_ax3 is not None else [lc_ax1, lc_ax2]

    if title is not None:
        fig.suptitle(title, x=0.5, y=1, va="bottom")

    return fig, axes


########################################
# Lightcurve optimization              #
########################################


def stack_fluxtable_rolling(
    fluxtable: Table,
    *,
    window_hours: float = 2.0,
) -> Table:
    """
    Rolling binning within each band: consecutive points are stacked as long as
    the time gap to the previous point is <= window_hours. It is used to produce
    a smoothed foreground lightcurve while keeping the original points available
    as a faint background layer.

    Stacking is performed in flux space using inverse-variance weights.
    time and (if present) zp are also weighted by the same weights.

    Returns a new Table with columns: time, flux, fluxerr, band, and zp if present.
    """
    if fluxtable is None or len(fluxtable) == 0:
        return fluxtable

    dt_max = float(window_hours) / 24.0  # hours -> days (JD)

    has_zp = "zp" in fluxtable.colnames
    out_rows: list[dict[str, Any]] = []

    # Ensure we don't modify the original table ordering/columns unexpectedly
    bands = np.unique(np.asarray(fluxtable["band"]))

    for band in bands:
        tab = fluxtable[fluxtable["band"] == band]
        if len(tab) == 0:
            continue

        t = np.asarray(tab["time"], dtype=float)
        f = np.asarray(tab["flux"], dtype=float)
        fe = np.asarray(tab["fluxerr"], dtype=float)

        # Sort by time
        idx = np.argsort(t)
        t, f, fe = t[idx], f[idx], fe[idx]

        zps = np.asarray(tab["zp"], dtype=float)[idx] if has_zp else None

        # Build bins of consecutive points with gap constraint
        start = 0
        n = len(t)
        while start < n:
            end = start + 1
            while end < n and (t[end] - t[end - 1]) <= dt_max:
                end += 1

            tt = t[start:end]
            ff = f[start:end]
            ee = fe[start:end]

            # Inverse-variance weights; guard against 0/NaN
            w = np.zeros_like(ee, dtype=float)
            ok = np.isfinite(ee) & (ee > 0) & np.isfinite(ff) & np.isfinite(tt)
            w[ok] = 1.0 / (ee[ok] ** 2)

            if np.sum(w) <= 0:
                # Fallback: unweighted mean if weights unusable
                t_mean = float(np.nanmean(tt))
                f_mean = float(np.nanmean(ff))
                # crude fallback error: standard error of mean if possible, else nan
                f_err = (
                    float(np.nanstd(ff) / np.sqrt(np.sum(np.isfinite(ff))))
                    if np.sum(np.isfinite(ff)) > 1
                    else float("nan")
                )
                row: dict[str, Any] = {
                    "time": t_mean,
                    "flux": f_mean,
                    "fluxerr": f_err,
                    "band": str(band),
                }
                if has_zp and zps is not None:
                    row["zp"] = float(np.nanmean(zps[start:end]))
                out_rows.append(row)
                start = end
                continue

            wsum = float(np.sum(w))
            t_mean = float(np.sum(w * tt) / wsum)
            f_mean = float(np.sum(w * ff) / wsum)
            f_err = float(np.sqrt(1.0 / wsum))

            row = {"time": t_mean, "flux": f_mean, "fluxerr": f_err, "band": str(band)}
            if has_zp and zps is not None:
                row["zp"] = float(np.sum(w * zps[start:end]) / wsum)

            out_rows.append(row)
            start = end

    # Keep dtype reasonably aligned
    cols: dict[str, Any] = {
        "time": [r["time"] for r in out_rows],
        "flux": [r["flux"] for r in out_rows],
        "fluxerr": [r["fluxerr"] for r in out_rows],
        "band": [r["band"] for r in out_rows],
    }
    if "zp" in (out_rows[0].keys() if out_rows else []):
        cols["zp"] = [r.get("zp", np.nan) for r in out_rows]

    return Table(cols)


def get_upperlimit_table(dps: Iterable[DataPoint] | None) -> Table | None:
    """
    Build a table of currently supported upper limits from datapoints.

    At present, this only supports ZTF upper limits with `diffmaglim` in the datapoint body.
    """
    if not dps:
        return None

    def filter_limits(dps_: Iterable[DataPoint]) -> list[dict]:
        return [
            dp["body"]
            for dp in dps_
            if "ZTF" in dp["tag"] and "diffmaglim" in dp["body"]
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

    return Table(
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


########################################
# Cutout image handling                #
########################################


def load_cutout_image(
    data_bytes: bytes,
    *,
    return_header: bool = False,
) -> np.ndarray | tuple[np.ndarray, fits.Header]:
    """
    Decode a cutout FITS image from bytes and return a 2D float array.

    Supports both:
      - ZTF cutouts: gzip-compressed FITS payload
      - LSST cutouts: plain FITS payload

    If return_header=True, also return the FITS header of the chosen 2D HDU.
    """
    # 1) Decompress if needed (ZTF), otherwise treat as plain FITS (LSST)
    try:
        with gzip.open(io.BytesIO(data_bytes), "rb") as f:
            payload = f.read()
    except BadGzipFile:
        payload = data_bytes

    with fits.open(io.BytesIO(payload), ignore_missing_simple=True) as hdul:
        # 2) Find first 2D image HDU
        data2d: np.ndarray | None = None
        hdr: fits.Header | None = None

        for hdu in hdul:
            if getattr(hdu, "data", None) is None:
                continue
            arr = np.squeeze(hdu.data)
            if isinstance(arr, np.ndarray) and arr.ndim == 2:
                data2d = arr.astype(float, copy=False)
                hdr = hdu.header
                break

        if data2d is None or hdr is None:
            raise ValueError("No 2D image HDU found in cutout FITS payload")

        if return_header:
            return data2d, hdr
        return data2d


def cutout_fov_arcsec(data2d: np.ndarray, hdr: fits.Header) -> float | None:
    """
    Calculate the Field of View (FoV) in arcseconds for a given cutout FITS image.
    """
    ny, nx = data2d.shape[-2], data2d.shape[-1]

    try:
        w = WCS(hdr)
        # pixel scales in deg/pix (for each axis)
        scales_deg = proj_plane_pixel_scales(w)  # ndarray, len>=2
        pixscale_arcsec = float(np.mean(scales_deg[:2]) * 3600.0)
        return pixscale_arcsec * float(max(nx, ny))
    except Exception:
        return None


def get_cached_cutout_image(
    data_bytes: bytes,
    *,
    cache_dir: str | None,
    cache_key: str,
    cutout_type: str,
) -> tuple[np.ndarray, float | None]:
    """
    Return a display-ready 2D cutout image and its FoV in arcsec.

    The cache stores the already processed image after:
      - FITS decoding
      - optional rotation
      - stretch / normalization for imshow
    """
    cache_path = None
    if cache_dir is not None:
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = os.path.join(cache_dir, f"{cache_key}_{cutout_type}.npz")
        if os.path.isfile(cache_path):
            try:
                cached = np.load(cache_path, allow_pickle=True)
                data = np.asarray(cached["data"], dtype=float)
                fov = cached["fov"]
                fov_val: float | None
                if np.ndim(fov) == 0:
                    fov_scalar = float(fov)
                    fov_val = fov_scalar if np.isfinite(fov_scalar) else None
                else:
                    fov_val = None
                return data, fov_val
            except Exception:
                pass

    data, hdr = load_cutout_image(data_bytes, return_header=True)
    cutout_fov = cutout_fov_arcsec(data, hdr)

    rotpa = hdr.get("ROTPA", None)
    if rotpa is not None:
        try:
            angle = -float(rotpa)
            data = ndimage.rotate(
                data,
                angle=angle,
                reshape=False,
                order=1,
                mode="constant",
                cval=np.nan,
                prefilter=False,
            )
        except Exception:
            pass

    finite = np.isfinite(data)
    if not np.any(finite):
        data_ = np.full_like(data, np.nan, dtype=float)
    else:
        vmin, vmax = np.percentile(data[finite], [50, 99.5])
        denom = vmax - vmin
        if not np.isfinite(denom) or denom <= 0:
            denom = 1.0
        data_ = visualization.AsinhStretch()((data - vmin) / denom)

    if cache_path is not None:
        with suppress(Exception):
            np.savez(
                cache_path,
                data=np.asarray(data_, dtype=float),
                fov=np.nan if cutout_fov is None else float(cutout_fov),
            )

    return np.asarray(data_, dtype=float), cutout_fov


def create_stamp_plot(
    cutouts: Mapping[str, Mapping[str, bytes]],
    ax,
    cutout_type: str,
    *,
    cache_dir: str | None = None,
    cache_key: str | None = None,
) -> float | None:
    """
    Render a Science/Template/Difference cutout into the provided axes and
    return an FOV in arcseconds to use as the finder FOV.
    """
    data_bytes = next(iter(cutouts.values()))[f"cutout{cutout_type}"]

    if cache_key is not None:
        data_, cutout_fov = get_cached_cutout_image(
            data_bytes,
            cache_dir=cache_dir,
            cache_key=cache_key,
            cutout_type=cutout_type,
        )
    else:
        data, hdr = load_cutout_image(data_bytes, return_header=True)
        cutout_fov = cutout_fov_arcsec(data, hdr)

        rotpa = hdr.get("ROTPA", None)
        if rotpa is not None:
            try:
                angle = -float(rotpa)
                data = ndimage.rotate(
                    data,
                    angle=angle,
                    reshape=False,
                    order=1,
                    mode="constant",
                    cval=np.nan,
                    prefilter=False,
                )
            except Exception:
                pass

        finite = np.isfinite(data)
        if not np.any(finite):
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(cutout_type, fontdict={"fontsize": "small"})
            return None

        vmin, vmax = np.percentile(data[finite], [50, 99.5])
        denom = vmax - vmin
        if not np.isfinite(denom) or denom <= 0:
            denom = 1.0
        data_ = visualization.AsinhStretch()((data - vmin) / denom)

    finite2 = np.isfinite(data_)
    if not np.any(finite2):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(cutout_type, fontdict={"fontsize": "small"})
        return cutout_fov

    ax.imshow(
        data_,
        norm=Normalize(*np.percentile(data_[finite2], [0.5, 99.5])),
        aspect="equal",
        cmap="viridis",
        origin="lower",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(cutout_type, fontdict={"fontsize": "small"})
    return cutout_fov


########################################
# Finder image handling                #
########################################


def render_finder_stamp(
    ax,
    ra: float,
    dec: float,
    *,
    cache_dir: str | None,
    name: str,
    size: int = 240,
    fov_arcsec: float = 45,
    crosshair_gap_frac: float = 0.1,
    matches: list[dict[str, Any]] | None = None,
    legend_ax=None,
) -> None:
    """
    Render a finder stamp.

    This includes:
    - Loading a HiPS finder stamp from a cache directory or downloading it.
    - Adding a gap crosshair to the finder stamp to indicate the center.
    - Grouping matches by their distance and plotting them as markers on the finder stamp.
    - Adding a legend with the grouped matches.
    """
    try:
        if cache_dir is not None:
            cache_key = f"{name}-{ra:.4f}-{dec:.4f}-{size}-{fov_arcsec:.4f}"
            cache_path = os.path.join(cache_dir, f"{cache_key}_FINDER.npz")
            if os.path.isfile(cache_path):
                cached = np.load(cache_path, allow_pickle=True)
                stamp = cached["data"]
                label = str(cached["label"])
            else:
                stamp, label = get_finder_stamp(
                    ra, dec, size=size, fov_arcsec=fov_arcsec
                )
                if stamp is not None and label is not None:
                    np.savez(
                        cache_path, data=np.asarray(stamp, dtype=float), label=label
                    )
        else:
            stamp, label = get_finder_stamp(ra, dec, size=size, fov_arcsec=fov_arcsec)
    except Exception:
        stamp, label = None, None

    if stamp is not None:
        half_fov = fov_arcsec / 2
        ax.imshow(
            stamp,
            aspect="equal",
            origin="upper",
            extent=[-half_fov, half_fov, -half_fov, half_fov],
        )
        add_gap_crosshair(ax, gap_frac=crosshair_gap_frac, lw=0.5, alpha=0.9)
        if matches:
            groups = group_matches(matches, ra, dec)
            for m in matches:
                if m["label"] == "T2NuclearFilter":
                    groups.append([m])

            plot_grouped_markers(
                ax,
                groups,
                center_ra=ra,
                center_dec=dec,
            )

            if legend_ax is not None:
                add_match_legend(legend_ax, groups)
        ax.set_title(label or "Finder", fontdict={"fontsize": "small"})
    else:
        ax.text(
            0.5, 0.5, "Finder unavailable", ha="center", va="center", fontsize="small"
        )
        ax.set_title("Finder", fontdict={"fontsize": "small"})

    ax.tick_params(axis="both", which="major", labelsize=6)
    ax.tick_params(axis="both", which="minor", labelsize=6)


def group_matches(matches, center_ra, center_dec, tol_arcsec=0.4):
    """
    Group matches that are closer than tol_arcsec on sky.
    Returns list of groups: each group is list of match dicts.
    """
    groups: list[list[dict[str, Any]]] = []

    cosdec = np.cos(np.deg2rad(center_dec))

    for m in matches:
        if m["label"] == "T2NuclearFilter":
            continue
        try:
            ra = float(m["ra"])
            dec = float(m["dec"])
        except Exception:
            continue

        placed = False

        for g in groups:
            ref = g[0]
            dra = (ra - ref["ra"]) * cosdec * 3600.0
            ddec = (dec - ref["dec"]) * 3600.0
            dist = np.hypot(dra, ddec)

            if dist < tol_arcsec:
                g.append(m)
                placed = True
                break

        if not placed:
            groups.append([m])

    return groups


COLORS = ["orangered", "lime", "magenta", "aqua", "hotpink"]

SHORT_MATCH_LABELS = {
    "GLADEv23": "GLADE",
    "NEDLVS": "NED LVS",
    "NEDz": "NED z",
    "NEDz_extcats": "NED ext",
    "SDSS_spec": "SDSS spec",
    "LSPhotoZZou": "LS Zou",
    "PS1_photoz": "PS1",
    "twoMPZ": "2MPZ",
    "wiseScosPhotoz": "WISE Scos",
    "wise_color": "WISE color",
    "milliquas": "Milliquas",
    "LSPhotoZTap": "LS DR10",
    "T2LSPhotoZTap": "LS DR10 Tap",
    "NED2026": "NED",
}

DIGEST_REDSHIFT_CATALOGS = [
    "NED2026",
    "SDSS_spec",
    "NEDz",
    "NEDz_extcats",
    "GLADEv23",
    "LSPhotoZZou",
    "twoMPZ",
    "PS1_photoz",
    "wiseScosPhotoz",
]


def short_match_label(label: str) -> str:
    return SHORT_MATCH_LABELS.get(label, label)


def plot_grouped_markers(
    ax,
    groups,
    *,
    center_ra,
    center_dec,
):
    """
    Project grouped RA/Dec positions into finder-image pixel coordinates and draw the
    colored circular markers. It also applies the RA-direction flip needed for the displayed
    finder orientation.
    """
    cosdec = np.cos(np.deg2rad(center_dec))

    for i, group in enumerate(groups):
        color = COLORS[i % len(COLORS)]

        # Use the first match as the representative position
        m = group[0]
        ra = float(m["ra"])
        dec = float(m["dec"])

        dra = (ra - center_ra) * cosdec * 3600
        ddec = (dec - center_dec) * 3600

        if not (np.isfinite(dra) and np.isfinite(ddec)):
            continue

        ax.plot(
            -dra,
            ddec,
            marker="o",
            markersize=6,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.5,
            zorder=5,
        )


def add_match_legend(ax, groups):
    """
    Draw a minimal text legend in a separate axis next to the finder. Each group gets one color.
    """
    ax.axis("off")

    y = 0.87
    line_step = 0.075
    block_gap = 0.035

    for i, group in enumerate(groups):
        color = COLORS[i % len(COLORS)]
        labels = sorted({short_match_label(str(m["label"])) for m in group})

        if not labels:
            continue

        # first line with color square
        ax.text(
            0.00,
            y,
            "■",
            color=color,
            fontsize=8,
            ha="left",
            va="top",
            family="monospace",
            alpha=0.95,
            transform=ax.transAxes,
        )
        ax.text(
            0.10,
            y,
            labels[0],
            color="black",
            fontsize=7,
            ha="left",
            va="top",
            family="monospace",
            alpha=0.9,
            transform=ax.transAxes,
        )
        y -= line_step

        # remaining labels underneath
        for lbl in labels[1:]:
            ax.text(
                0.10,
                y,
                lbl,
                color="black",
                fontsize=7,
                ha="left",
                va="top",
                family="monospace",
                alpha=0.9,
                transform=ax.transAxes,
            )
            y -= line_step

        y -= block_gap

        if y < 0.05:
            break


def get_finder_stamp(
    ra: float,
    dec: float,
    size: int = 240,
    surveys: list[str] | None = None,
    timeout: float = 8.0,
    fov_arcsec: float = 45.0,
) -> tuple[np.ndarray | None, str | None]:
    """
    Download a finder stamp centered at the given RA and Dec.

    The stamp is fetched from CDS HiPS services as a JPEG image and returned
    as a 2D grayscale numpy array together with a survey label.
    """

    if surveys is None:
        surveys = [
            "CDS/P/DESI-Legacy-Surveys/DR10/color",
            "CDS/P/PanSTARRS/DR1/color-z-zg-g",
            "CDS/P/DECaLS/DR5/color",
            "CDS/P/DES-DR2/ColorIRG",
            "CDS/P/Skymapper/DR4/color",
            "CDS/P/DSS2/color",
        ]

    labels = {
        "CDS/P/DESI-Legacy-Surveys/DR10/color": "LS DR10",
        "CDS/P/DES-DR2/ColorIRG": "DES",
        "CDS/P/Skymapper/DR4/color": "SkyMapper",
        "CDS/P/DECaLS/DR5/color": "DECaLS",
        "CDS/P/DSS2/color": "DSS2",
        "CDS/P/PanSTARRS/DR1/color-z-zg-g": "PS1",
    }

    BASE_URLS = [
        "https://alasky.cds.unistra.fr/hips-image-services/hips2fits",
        "https://alaskybis.cds.unistra.fr/hips-image-services/hips2fits",
    ]

    fov_deg = float(fov_arcsec) / 3600.0

    try:
        from PIL import Image, ImageOps  # noqa: PLC0415
    except Exception:
        return None, None

    for hips in surveys:
        try:
            params: dict[str, str | int | float] = {
                "hips": hips,
                "ra": ra,
                "dec": dec,
                "fov": fov_deg,
                "width": size,
                "height": size,
                "format": "jpg",
                "projection": "TAN",
            }
            r = None
            for base_url in BASE_URLS:
                try:
                    r = requests.get(base_url, params=params, timeout=timeout)
                except Exception:
                    continue

                if r.status_code == 200 and r.content:
                    break

                r = None

            if r is None:
                continue

            img = Image.open(io.BytesIO(r.content))
            img = ImageOps.exif_transpose(img).convert("RGB")
            arr = np.asarray(img).astype(float) / 255.0

            return arr, labels.get(hips, hips)

        except Exception:
            continue

    return None, None


def add_gap_crosshair(
    ax, *, gap_frac: float = 0.2, lw: float = 1.0, alpha: float = 0.8
) -> None:
    """
    Draw a crosshair in axes coordinates that leaves a central gap.
    gap_frac=0.2 leaves the central 20% of the axis free in both directions.
    """
    if not (0.0 < gap_frac < 1.0):
        return

    g0 = 0.5 - gap_frac / 2.0
    g1 = 0.5 + gap_frac / 2.0

    # Vertical segments (x=0.5), leaving a gap [g0, g1]
    ax.plot(
        [0.5, 0.5], [0.0, g0], transform=ax.transAxes, lw=lw, alpha=alpha, color="white"
    )
    ax.plot(
        [0.5, 0.5], [g1, 1.0], transform=ax.transAxes, lw=lw, alpha=alpha, color="white"
    )

    # Horizontal segments (y=0.5)
    ax.plot(
        [0.0, g0], [0.5, 0.5], transform=ax.transAxes, lw=lw, alpha=alpha, color="white"
    )
    ax.plot(
        [g1, 1.0], [0.5, 0.5], transform=ax.transAxes, lw=lw, alpha=alpha, color="white"
    )


class PlotNuclearFilterLightcurves(AbsPhotoT3Unit, AbsTabulatedT2Unit):
    """
    Create a (pdf) plot summarizing lightcurves of candidates provided to the unit.

    Features:
      - Optional stacked datapoints
      - Optional model lightcurves
      - Optional thumbnails (ZTFCutoutImages / LSSTCutoutImages)
      - Optional finder stamp including markers of closest match
      - Optional link to TNS (from TNSReports complement)
      - Optional Fritz link (disabled for LSST-only sources)
    """

    pdf_path: None | str = None  # Will create random if not set
    titleprefix: str = "AMPEL: "
    save_png: bool = False
    finder_cache_dir: str | None = "./finder_cache"
    cutout_cache_dir: str | None = "./cutout_cache"
    save_dir: str = "./images"

    include_cutouts: bool = False
    internal_cutout_fetch: bool = False
    cutout_eligible: str = "brightest"  # first | last | brightest
    cutout_fallback_tries: int = 2
    cutout_fallback_min_separation_hours: float = 15.0
    lsst_archive_insecure: bool = True
    lsst_archive_url: str = "https://ampel-dev.ia.zeuthen.desy.de/api/lsst/archive/v1/"
    ztf_archive_url: str | None = None

    fritzlink: bool = True

    stacking: bool = False
    stacking_window_hours: float = 2.0
    stacking_alpha: float = 0.22

    local_catalogs: list[str] = []  # Names of local catalogs to include

    def post_init(self) -> None:
        """
        Post-initialization routine for PlotTransientLightcurves.

        Creates save directories (including image cache if requested), sets up
        Slack WebClient (if enabled), and initializes ParsnipSncosmoSource
        if include_model_lightcurves is True.
        """
        os.makedirs(self.save_dir, exist_ok=True)
        if self.finder_cache_dir:
            os.makedirs(self.finder_cache_dir, exist_ok=True)
        if self.cutout_cache_dir:
            os.makedirs(self.cutout_cache_dir, exist_ok=True)
        if not self.pdf_path:
            self.pdf_path = tempfile.mkstemp(".pdf", "candidates", self.save_dir)[1]

        if self.ztf_archive_url is None:
            try:
                context = getattr(self, "context", None)
                if context is not None:
                    self.ztf_archive_url = context.config.get(
                        "resource.ampel-ztf/archive",
                        str,
                        raise_exc=True,
                    )
                else:
                    self.ztf_archive_url = None
            except Exception:
                self.ztf_archive_url = None

        self._http_session = requests.Session()

        self._ztf_cutout_session = None
        if self.ztf_archive_url is not None:
            self._ztf_cutout_session = BaseUrlSession(base_url=self.ztf_archive_url)

    ########################################
    # Collecting information from T2       #
    ########################################

    def photz_from_t2s(self, tview: TransientView) -> list:
        """
        Collect information from T2CatalogMatch, T2LSPhotoZ and T2CatalogMatchLocal documents, return as list of str.
        """
        photz_list = []
        t2res = tview.get_t2_body(unit="T2DigestRedshifts")

        if isinstance(t2res, dict) and t2res.get("ampel_z", -10) > 0:
            t2cat = tview.get_t2_body(unit="T2CatalogMatch")
            if isinstance(t2cat, dict):
                if t2cat.get("GLADEv23") and t2cat["GLADEv23"].get("z") is not None:
                    photz_list.append("GLADE: {:.2f}".format(t2cat["GLADEv23"]["z"]))

                if t2cat.get("SDSS_spec") and t2cat["SDSS_spec"].get("z") is not None:
                    photz_list.append("SDSS: {:.2f}".format(t2cat["SDSS_spec"]["z"]))

                if t2cat.get("NED2026") and t2cat["NED2026"].get("z") is not None:
                    z = float(t2cat["NED2026"]["z"])
                    zerr = t2cat["NED2026"].get("zunc")
                    if zerr is not None:
                        photz_list.append(f"NED: {z:.4f} ± {float(zerr):.4f}")
                    else:
                        photz_list.append(f"NED: {z:.4f}")

                if t2cat.get("NEDz") and t2cat["NEDz"].get("z") is not None:
                    z = float(t2cat["NEDz"]["z"])
                    zerr = t2cat["NEDz"].get("zunc")
                    if zerr is not None:
                        photz_list.append(f"NED z: {z:.4f} ± {float(zerr):.4f}")
                    else:
                        photz_list.append(f"NED z: {z:.4f}")

                if (
                    t2cat.get("NEDz_extcats")
                    and t2cat["NEDz_extcats"].get("z") is not None
                ):
                    z = float(t2cat["NEDz_extcats"]["z"])
                    zerr = t2cat["NEDz_extcats"].get("zunc")
                    if zerr is not None:
                        photz_list.append(f"NED ext: {z:.4f} ± {float(zerr):.4f}")
                    else:
                        photz_list.append(f"NED ext: {z:.4f}")

                if (
                    t2cat.get("wiseScosPhotoz")
                    and t2cat["wiseScosPhotoz"].get("zPhoto_Corr") is not None
                ):
                    photz_list.append(
                        "wiseScos: {:.2f}".format(
                            float(t2cat["wiseScosPhotoz"]["zPhoto_Corr"])
                        )
                    )

                if (
                    t2cat.get("PS1_photoz")
                    and t2cat["PS1_photoz"].get("z_phot") is not None
                ):
                    ps1_photz = float(t2cat["PS1_photoz"]["z_phot"]) / 1000
                    zerr_raw = t2cat["PS1_photoz"].get("z_photErr")
                    if zerr_raw is not None:
                        ps1_photzerr = float(zerr_raw) / 10000
                        photz_list.append(f"PS1: {ps1_photz:.2f} ± {ps1_photzerr:.2f}")
                    else:
                        photz_list.append(f"PS1: {ps1_photz:.2f}")

                if t2cat.get("LSPhotoZZou"):
                    if (
                        t2cat["LSPhotoZZou"].get("specz") is not None
                        and t2cat["LSPhotoZZou"]["specz"] > -0.1
                    ):
                        photz_list.append(
                            f"LSZou spec: {float(t2cat['LSPhotoZZou']['specz']):.4f}"
                        )

                    if t2cat["LSPhotoZZou"].get("photoz") is not None:
                        z = float(t2cat["LSPhotoZZou"]["photoz"])
                        zerr = t2cat["LSPhotoZZou"].get("e_photoz")
                        if zerr is not None:
                            photz_list.append(f"LSZou: {z:.2f} ± {float(zerr):.2f}")
                        else:
                            photz_list.append(f"LSZou: {z:.2f}")

                if t2cat.get("twoMPZ"):
                    if (
                        t2cat["twoMPZ"].get("zSpec") is not None
                        and t2cat["twoMPZ"]["zSpec"] > -0.1
                    ):
                        photz_list.append(
                            f"2MPZ spec: {float(t2cat['twoMPZ']['zSpec']):.4f}"
                        )

                    if (
                        t2cat["twoMPZ"].get("zPhoto") is not None
                        and t2cat["twoMPZ"]["zPhoto"] > -0.1
                    ):
                        photz_list.append(
                            f"2MPZ: {float(t2cat['twoMPZ']['zPhoto']):.2f}"
                        )

            t2ls = tview.get_t2_body(unit="T2LSPhotoZTap")
            if isinstance(t2ls, dict) and t2ls.get("T2LSPhotoZTap"):
                if (
                    t2ls["T2LSPhotoZTap"].get("z_spec") is not None
                    and t2ls["T2LSPhotoZTap"]["z_spec"] > -1
                ):
                    photz_list.append(
                        f"LS DR10 spec: {float(t2ls['T2LSPhotoZTap']['z_spec']):.4f}"
                    )

                if t2ls["T2LSPhotoZTap"].get("z_phot_median") is not None:
                    z = float(t2ls["T2LSPhotoZTap"]["z_phot_median"])
                    zerr = t2ls["T2LSPhotoZTap"].get("z_phot_std")
                    if zerr is not None:
                        photz_list.append(f"LS DR10: {z:.2f} ± {float(zerr):.2f}")
                    else:
                        photz_list.append(f"LS DR10: {z:.2f}")

            t2cat_local = tview.get_t2_body(unit="T2CatalogMatchLocal")
            if isinstance(t2cat_local, dict):
                for catalog_name in self.local_catalogs:
                    entry = t2cat_local.get(catalog_name)
                    if not isinstance(entry, Mapping):
                        continue
                    if entry.get("z") is None:
                        continue

                    z = float(entry["z"])
                    zerr = (
                        entry.get("zunc")
                        or entry.get("z_err")
                        or entry.get("zErr")
                        or entry.get("z_photErr")
                        or entry.get("e_photoz")
                        or entry.get("zPhotoErr")
                    )
                    if zerr is not None:
                        photz_list.append(
                            f"{catalog_name} (local): {z:.4f} ± {float(zerr):.4f}"
                        )
                    else:
                        photz_list.append(f"{catalog_name} (local): {z:.4f}")

        return photz_list

    def attributes_from_t2(
        self,
        tview: TransientView,
        nearby_z: float = 0.02,
    ) -> tuple[list, Any]:
        """
        Collect information from potential T2 documents, return as list of str.
        """
        attributes: list[str] = []
        z = None

        t2res = tview.get_t2_body(unit="T2DigestRedshifts")
        if isinstance(t2res, dict) and t2res.get("ampel_z", -10) > 0:
            attributes.append("Ampel z: {:.2f}".format(t2res["ampel_z"]))
            z = t2res["ampel_z"]
            if t2res.get("ampel_z", 999) < nearby_z:
                attributes.append("Nearby")
            dist2host = t2res["ampel_dist"]
            attributes.append(f"Dist2host: {dist2host:.2f} arcsec")

        return attributes, z

    def finder_matches_from_t2(self, tview: TransientView) -> list[dict[str, Any]]:
        """
        Collect sky positions of catalog matches that can be overplotted in the finder stamp.
        """
        out: list[dict[str, Any]] = []

        # catalogs that are most relevant as possible host / counterpart positions
        # Keep this in sync with T2DigestRedshifts._get_catalogmatch_groupz()
        preferred_catalogs = DIGEST_REDSHIFT_CATALOGS

        unit_mapping = get_catalog_position_unit_map()

        # minimal entry to mock LS data structure results from T2LSPhotoZTap
        unit_mapping["T2LSPhotoZTap"] = {"ra": {"unit": "deg"}}

        def _append_match(cat_name: str, entry: Mapping[str, Any]) -> None:
            ra_key = None
            dec_key = None

            for rk in ("ra", "RA", "raMean"):
                if rk in entry and entry[rk] is not None:
                    ra_key = rk
                    break

            for dk in ("dec", "Dec", "decMean"):
                if dk in entry and entry[dk] is not None:
                    dec_key = dk
                    break

            if ra_key is None or dec_key is None:
                return

            try:
                ra_val = float(entry[ra_key])
                dec_val = float(entry[dec_key])
            except Exception:
                return

            if not np.isfinite(ra_val) or not np.isfinite(dec_val):
                return

            label = cat_name
            if cat_name == "T2LSPhotoZTap":
                label = "LSPhotoZTap"

            if (u := unit_mapping[cat_name]["ra"]["unit"]) is not None and (
                u.lower().startswith("rad")
            ):
                ra_val = np.degrees(ra_val)
                dec_val = np.degrees(dec_val)

            item: dict[str, Any] = {
                "ra": ra_val,
                "dec": dec_val,
                "label": label,
            }

            if "dist2transient" in entry:
                with suppress(Exception):
                    item["dist2transient"] = float(entry["dist2transient"])

            out.append(item)

        t2cat = tview.get_t2_body(unit="T2CatalogMatch")
        if isinstance(t2cat, dict):
            for cat_name in preferred_catalogs:
                entry = t2cat.get(cat_name)
                if isinstance(entry, Mapping):
                    _append_match(cat_name, entry)

        t2cat_local = tview.get_t2_body(unit="T2CatalogMatchLocal")
        if isinstance(t2cat_local, dict):
            for cat_name in self.local_catalogs:
                entry = t2cat_local.get(cat_name)
                if isinstance(entry, Mapping):
                    _append_match(cat_name, entry)

        t2ls = tview.get_t2_body(unit="T2LSPhotoZTap")
        if isinstance(t2ls, dict) and isinstance(t2ls.get("T2LSPhotoZTap"), Mapping):
            _append_match("T2LSPhotoZTap", t2ls["T2LSPhotoZTap"])

        # remove near-duplicates by sky position + label
        deduped: list[dict[str, Any]] = []
        seen: set[tuple[str, int, int]] = set()

        for item in out:
            key = (
                str(item["label"]),
                round(float(item["ra"]) * 1e6),
                round(float(item["dec"]) * 1e6),
            )
            if key in seen:
                continue
            seen.add(key)
            deduped.append(item)

        return deduped

    ########################################
    # Collecting photometry                #
    ########################################

    # Because TransientView does not guarantee any id convention, we introduce these helpers
    # to treat all datapoints as possible photopoints/upperlimits.
    # Eventually, we need to fix that with a proper tabulator-like unit!

    def _iter_t0(self, tview: TransientView) -> Sequence[DataPoint]:
        return tview.t0 or []

    def _is_flux_dp(self, dp: DataPoint) -> bool:
        """
        Checks if a given datapoint can be considered a flux-like measurement.
        """
        body = dp.get("body") or {}
        if not isinstance(body, dict):
            return False

        # Common flux-like measurement fields
        for fmeas in ("flux", "fluxpsf", "psfFlux", "mag", "magpsf"):
            if body.get(fmeas) is not None:
                return True

        return False

    def _is_supported_upperlimit_dp(self, dp: DataPoint) -> bool:
        """
        Checks if a given datapoint can be considered a supported upper limit.
        Currently, only ZTF-style upper limits are supported by get_upperlimit_table().
        """
        body = dp.get("body") or {}
        tags = dp.get("tag") or []

        if not isinstance(body, dict):
            return False

        if isinstance(tags, str):
            tags = [tags]

        # Currently only ZTF-style upper limits are supported by get_upperlimit_table()
        return "diffmaglim" in body and "ZTF" in tags

    def _get_photopoints_any_id(self, tview: TransientView) -> list[DataPoint]:
        dps = self._iter_t0(tview)
        return [dp for dp in dps if self._is_flux_dp(dp)]

    def _get_upperlimits_any_id(self, tview: TransientView) -> list[DataPoint]:
        dps = self._iter_t0(tview)
        return [dp for dp in dps if self._is_supported_upperlimit_dp(dp)]

    ########################################
    # Retrieving cutouts internally        #
    ########################################

    # To speed up local processing, we need to cache the cutouts. For that, we need to incorporate the funcionalities of
    # LSSTCutoutImages and ZTFCutoutImages. The original units are still usable as complements.

    def _get_raw_cutout_cache_path(self, survey: str, candidate_id: int) -> str | None:
        """
        Returns the path to the cache directory for raw cutouts associated with a given survey and candidate ID.
        """
        if not self.cutout_cache_dir:
            return None

        raw_dir = os.path.join(self.cutout_cache_dir, "raw")
        os.makedirs(raw_dir, exist_ok=True)
        return os.path.join(raw_dir, f"{survey}_{candidate_id}.npz")

    def _select_lsst_cutout_candidates(
        self,
        photopoints: Sequence[DataPoint],
    ) -> list[int]:
        """
        Selects the LSST cutout candidates from a list of photopoints.

        If cutout_eligible == "last", the candidates are returned in reverse order.
        If cutout_eligible == "first", the candidates are returned in their original order.
        If cutout_eligible == "brightest", the candidates are sorted by their psfFlux in descending order.
        Otherwise, the candidates are returned in reverse order.

        Candidates with a midpointMjdTai closer than cutout_fallback_min_separation_hours to the previously
        accepted candidates are skipped. In this way, the same observation run is not selected twice.
        """
        pps = sorted(
            [
                pp
                for pp in photopoints
                if "LSST_DP" in pp.get("tag", [])
                and isinstance(pp.get("body"), dict)
                and pp["body"].get("diaSourceId") is not None
                and pp["body"].get("midpointMjdTai") is not None
            ],
            key=lambda pp: pp["body"]["midpointMjdTai"],
        )

        if not pps:
            return []

        def _diasource_id(pp: DataPoint) -> int:
            body = pp.get("body", {})
            dsid = body.get("diaSourceId")
            if dsid is None:
                raise KeyError(
                    f"No diaSourceId in photopoint body keys={list(body.keys())}"
                )
            return int(dsid)

        if self.cutout_eligible == "last":
            ordered = list(reversed(pps))
        elif self.cutout_eligible == "first":
            ordered = pps
        elif self.cutout_eligible == "brightest":
            ordered = sorted(
                pps,
                key=lambda pp: float(pp["body"].get("psfFlux", float("-inf"))),
                reverse=True,
            )
        else:
            ordered = list(reversed(pps))

        sep_days = max(float(self.cutout_fallback_min_separation_hours), 0.0) / 24.0

        candids: list[int] = []
        accepted_times: list[float] = []
        seen_ids: set[int] = set()

        for pp in ordered:
            body = pp.get("body", {})

            try:
                candid = _diasource_id(pp)
                t = float(body["midpointMjdTai"])
            except Exception:
                continue

            if candid in seen_ids:
                continue

            if (
                accepted_times
                and sep_days > 0
                and min(abs(t - t0) for t0 in accepted_times) < sep_days
            ):
                continue

            seen_ids.add(candid)
            accepted_times.append(t)
            candids.append(candid)

        return candids

    def _select_ztf_cutout_candidates(
        self,
        photopoints: Sequence[DataPoint],
    ) -> list[int]:
        """
        Selects the ZTF cutout candidates. Analogous to _select_lsst_cutout_candidates.
        """
        pps = sorted(
            [
                pp
                for pp in photopoints
                if pp.get("id") is not None
                and isinstance(pp.get("body"), dict)
                and pp["id"] > 0
                and pp["body"].get("jd") is not None
            ],
            key=lambda pp: pp["body"]["jd"],
        )

        if not pps:
            return []

        if self.cutout_eligible == "last":
            candids = [int(pps[-1]["id"])]
        elif self.cutout_eligible == "first":
            candids = [int(pps[0]["id"])]
        elif self.cutout_eligible == "brightest":
            candids = [int(min(pps, key=lambda pp: pp["body"]["magpsf"])["id"])]
        else:  # all
            candids = [int(pp["id"]) for pp in pps]

        out: list[int] = []
        seen: set[int] = set()

        for candid in candids:
            if candid in seen:
                continue
            seen.add(candid)
            out.append(candid)

        return out

    @backoff.on_exception(
        backoff.expo,
        (requests.ConnectionError, requests.Timeout),
        max_tries=5,
        factor=2,
    )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: (
            not isinstance(e, requests.HTTPError)
            or e.response is None
            or e.response.status_code not in {502, 503, 504, 429, 408}
        ),
        max_time=60,
    )
    def _download_lsst_cutout(self, dia_source_id: int) -> dict[str, bytes] | None:
        """
        Download a cutout from the LSST archive.

        Behavior:
        - 404 -> no cutout available -> return None
        - temporary network / server problems -> retry
        - persistent timeout / server failure -> raise and stop job
        """
        response = self._http_session.get(
            f"{self.lsst_archive_url.rstrip('/')}/alert/{dia_source_id}/cutouts",
            verify=not self.lsst_archive_insecure,
            timeout=20,
        )

        if response.status_code == 404:
            return None

        response.raise_for_status()
        json = response.json()
        return {
            k: base64.b64decode(json[k])
            for k in ["cutoutScience", "cutoutTemplate", "cutoutDifference"]
        }

    @backoff.on_exception(
        backoff.expo,
        requests.ConnectionError,
        max_tries=5,
        factor=10,
    )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: (
            not isinstance(e, requests.HTTPError)
            or e.response is None
            or e.response.status_code not in {502, 503, 504, 429, 408}
        ),
        max_time=60,
    )
    def _download_ztf_cutout(self, candid: int) -> dict[str, bytes] | None:
        """
        Download a cutout from the ZTF archive. Uses BaseURLSession like ZTFCutoutImages.py.
        """
        if self._ztf_cutout_session is None:
            return None

        response = self._ztf_cutout_session.get(f"alert/{candid}/cutouts")

        if response.status_code == 404:
            return None

        response.raise_for_status()
        json = response.json()
        return {
            k: base64.b64decode(json[k]["stampData"])
            for k in ["cutoutScience", "cutoutTemplate", "cutoutDifference"]
        }

    def _get_lsst_cutout(self, dia_source_id: int) -> dict[str, bytes] | None:
        """
        Retrieve a cutout from the LSST archive, from the cache if possible
        """
        cache_path = self._get_raw_cutout_cache_path("LSST", dia_source_id)

        if cache_path and os.path.isfile(cache_path):
            try:
                cached = np.load(cache_path, allow_pickle=True)
                status = str(cached.get("status", "ok"))
                if status != "ok":
                    return None

                return {
                    "cutoutScience": bytes(cached["cutoutScience"].tolist()),
                    "cutoutTemplate": bytes(cached["cutoutTemplate"].tolist()),
                    "cutoutDifference": bytes(cached["cutoutDifference"].tolist()),
                }
            except Exception:
                pass

        cutout = self._download_lsst_cutout(dia_source_id)

        if cutout is None:
            if cache_path:
                with suppress(Exception):
                    np.savez(cache_path, status="missing")
            return None

        if cache_path:
            with suppress(Exception):
                np.savez(
                    cache_path,
                    status="ok",
                    cutoutScience=np.frombuffer(
                        cutout["cutoutScience"], dtype=np.uint8
                    ),
                    cutoutTemplate=np.frombuffer(
                        cutout["cutoutTemplate"], dtype=np.uint8
                    ),
                    cutoutDifference=np.frombuffer(
                        cutout["cutoutDifference"], dtype=np.uint8
                    ),
                )

        return cutout

    def _get_ztf_cutout(self, candid: int) -> dict[str, bytes] | None:
        """
        Retrieve a cutout from the ZTF archive. Analogous to _get_lsst_cutout.
        """
        cache_path = self._get_raw_cutout_cache_path("ZTF", candid)

        if cache_path and os.path.isfile(cache_path):
            try:
                cached = np.load(cache_path, allow_pickle=True)
                status = str(cached.get("status", "ok"))
                if status != "ok":
                    return None

                return {
                    "cutoutScience": bytes(cached["cutoutScience"].tolist()),
                    "cutoutTemplate": bytes(cached["cutoutTemplate"].tolist()),
                    "cutoutDifference": bytes(cached["cutoutDifference"].tolist()),
                }
            except Exception:
                pass

        cutout = self._download_ztf_cutout(candid)

        if cutout is None:
            if cache_path:
                with suppress(Exception):
                    np.savez(cache_path, status="missing")
            return None

        if cache_path:
            with suppress(Exception):
                np.savez(
                    cache_path,
                    status="ok",
                    cutoutScience=np.frombuffer(
                        cutout["cutoutScience"], dtype=np.uint8
                    ),
                    cutoutTemplate=np.frombuffer(
                        cutout["cutoutTemplate"], dtype=np.uint8
                    ),
                    cutoutDifference=np.frombuffer(
                        cutout["cutoutDifference"], dtype=np.uint8
                    ),
                )

        return cutout

    def _build_cutouts_from_photopoints(
        self,
        photopoints: Sequence[DataPoint],
    ) -> tuple[Mapping[str, Mapping[str, bytes]] | None, str | None]:
        """
        Build cutouts from a list of photopoints.

        This function tries to find a matching cutout for each photopoint in the list.
        It will return the first matching cutout it finds.
        """
        max_tries = max(int(self.cutout_fallback_tries), 1)

        lsst_candids = self._select_lsst_cutout_candidates(photopoints)
        for candid in lsst_candids[:max_tries]:
            cutout = self._get_lsst_cutout(candid)
            if cutout is not None:
                cache_key = f"LSST_{candid}"
                return {cache_key: cutout}, cache_key

        ztf_candids = self._select_ztf_cutout_candidates(photopoints)
        for candid in ztf_candids[:max_tries]:
            cutout = self._get_ztf_cutout(candid)
            if cutout is not None:
                cache_key = f"ZTF_{candid}"
                return {cache_key: cutout}, cache_key

        return None, None

    ########################################
    # The actual process of plotting       #
    ########################################

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None
    ) -> UBson | UnitResult:
        """
        This function loops through all transients, collects their data and parameters by using the
        prior methods, and feeds them to fig_from_fluxtable. It then saves the resulting figure as a PDF.
        Optionally, it uploads the PDF to Slack.
        """
        buffer = io.BytesIO()
        with PdfPages(self.pdf_path or buffer) as pdf:
            for tran_view in gen:
                photopoints = self._get_photopoints_any_id(tran_view)
                if not photopoints:
                    self.logger.debug("No photopoints", extra={"stock": tran_view.id})
                    continue

                sncosmo_table = self.get_flux_table(photopoints)
                if sncosmo_table is None or len(sncosmo_table) == 0:
                    self.logger.debug("No flux table", extra={"stock": tran_view.id})
                    continue

                sncosmo_table = sncosmo_table[sncosmo_table["flux"] > 0]
                if len(sncosmo_table) == 0:
                    self.logger.debug("No pos flux", extra={"stock": tran_view.id})
                    continue

                (ra, dec) = self.get_pos(photopoints)
                name = " ".join(map(str, self.get_stock_name(photopoints)))
                ampelid = tran_view.id

                upperlimits = self._get_upperlimits_any_id(tran_view)
                ulim_table = get_upperlimit_table(upperlimits)

                title = f"{self.titleprefix}: {name!r}-{ampelid!r} @ RA {ra:.3f} Dec {dec:.3f}"

                (attributes, z) = self.attributes_from_t2(tran_view)
                photz_list = self.photz_from_t2s(tran_view)

                tnsname = None
                if (
                    isinstance(tran_view.extra, dict)
                    and "TNSReports" in tran_view.extra
                ):
                    tnsname = next(iter(tran_view.extra["TNSReports"])).get(
                        "objname", None
                    )

                cutouts = None
                cutout_cache_key = None

                if self.include_cutouts:
                    if self.internal_cutout_fetch:
                        cutouts, cutout_cache_key = (
                            self._build_cutouts_from_photopoints(photopoints)
                        )
                        if cutouts:
                            self.logger.debug(
                                "Using internally fetched cutouts",
                                extra={
                                    "stock": tran_view.id,
                                    "cache_key": cutout_cache_key,
                                },
                            )
                    elif isinstance(tran_view.extra, dict):
                        for cutout_source in ("LSSTCutoutImages", "ZTFCutoutImages"):
                            raw = tran_view.extra.get(cutout_source)
                            if not isinstance(raw, dict):
                                continue
                            filtered = {k: v for k, v in raw.items() if v}
                            if filtered:
                                cutouts = filtered
                                try:
                                    cutout_cache_key = str(next(iter(filtered.keys())))
                                except Exception:
                                    cutout_cache_key = name
                                self.logger.debug(
                                    "Found complemented cutouts",
                                    extra={
                                        "stock": tran_view.id,
                                        "source": cutout_source,
                                    },
                                )
                                break

                finder_matches = self.finder_matches_from_t2(tran_view)

                fig, _axes = fig_from_fluxtable(
                    name,
                    str(ampelid),
                    ra,
                    dec,
                    sncosmo_table,
                    ulim_table,
                    nuclear_filter_res=tran_view.get_t2_body(unit="T2NuclearFilter"),
                    attributes=attributes,
                    photz_list=photz_list,
                    fritzlink=self.fritzlink,
                    tnsname=tnsname,
                    z=z,
                    cutouts=cutouts,
                    cutout_cache_key=cutout_cache_key,
                    finder_cache_dir=self.finder_cache_dir,
                    cutout_cache_dir=self.cutout_cache_dir,
                    stacking=self.stacking,
                    stacking_window_hours=self.stacking_window_hours,
                    stacking_alpha=self.stacking_alpha,
                    finder_matches=finder_matches,
                    title=title,
                )
                pdf.savefig(fig)
                if self.save_png:
                    plt.savefig(os.path.join(self.save_dir, str(tran_view.id) + ".png"))
                plt.close(fig)

        return None
