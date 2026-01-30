#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/PlotTransientLightcurves.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                11.06.2018
# Last Modified Date:  21.01.2026
# Last Modified By:    Felix Fischer <firstname.martin.lastname@desy.de>

import gzip
import io
import os
import tempfile
from collections.abc import Generator, Iterable, Mapping, Sequence
from gzip import BadGzipFile
from typing import Any

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
from scipy import ndimage
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
from ztfquery.utils.stamps import get_ps_stamp

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.secret.NamedSecret import NamedSecret
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
    title: None | str = None,
    tnsname: str | None = None,
    fritzlink: bool = True,
    attributes: list[str] | None = None,
    photz_list: list[str] | None = None,
    classprobs: Mapping[str, float] | None = None,
    cutouts: Mapping[str, Mapping[str, bytes]] | None = None,
    mag_range: None | list = None,
    z: float | None = None,
    zp: float = 25.0,
    legend: bool = True,
    grid_interval: None | int = None,
    t_0_jd: None | float = None,
    cutout_cache_dir: None | str = ".",
):
    """
    Create a lightcurve figure (in mag space) based on flux tables from AMPEL.

    If cutouts are provided, a multi-panel layout is used:
      - Science/Template/Difference cutouts
      - Finder stamp (PS1 preferred; HiPS fallback)
      - Lightcurve
      - Optional classification radar + text table
      - General info box (coords, attributes, photo-z)

    If no cutouts are provided, only the lightcurve is shown.
    """
    attributes = attributes or []
    photz_list = photz_list or []
    classprobs = dict(classprobs or {})

    cls_threshold = 0.01  # minimum class probability to be shown
    cls_items = filtered_classprobs(classprobs, threshold=cls_threshold)
    has_cutouts = bool(cutouts)
    has_classprobs = bool(cls_items)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Prefer per-point zp from tabulator output, fallback to function argument
    if "zp" in fluxtable.colnames:
        zps = np.asarray(fluxtable["zp"], dtype=float)
    else:
        zps = np.full(len(fluxtable), float(zp), dtype=float)

    # Transform to mags (AB)
    fluxtable["mag"] = -2.5 * np.log10(fluxtable["flux"]) + zps
    fluxtable["magerr"] = np.abs(
        -2.5 * fluxtable["fluxerr"] / (fluxtable["flux"] * np.log(10))
    )

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

    info.append("------------------------")
    if attributes:
        info.extend(attributes)
    if photz_list:
        info.append("------------------------")
        info.append("Photo-z:")
        info.extend(photz_list)

    fig = plt.figure(figsize=(8, 5))
    lc_ax3 = None

    # Helper functions for layout
    def _plot_cutouts_and_get_fov(*, cutouts, axes) -> float | None:
        """
        Render Science/Template/Difference cutouts into the provided axes and
        return an FOV in arcseconds to use as the finder FOV.
        """
        cutout_fov = None
        for ax_, typ in zip(axes, ["Science", "Template", "Difference"], strict=False):
            fov = create_stamp_plot(cutouts=cutouts, ax=ax_, cutout_type=typ)
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

    def _make_cls_lines(cls_items, max_items: int = 8) -> list[str]:
        """
        Format a compact text table of class probabilities for the attributes box.
        Only the top entries are shown.
        """
        lines = ["Class Probabilities (p>0.01):"]
        for k, v in cls_items[:max_items]:
            lines.append(f"{k}: {v:.3f}")
        return lines

    def _setup_gs_2x32(fig):
        """
        Create and space the common 2x32 GridSpec used for layouts with cutous.
        """
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
            top=0.93,
            bottom=0.10,
            wspace=0.35,
            hspace=0.45,
        )
        return gs

    # Layout A: Cutouts + finder stamp + classification radar.
    if has_cutouts and has_classprobs:
        assert cutouts is not None

        gs = _setup_gs_2x32(fig)

        # Top row: 3 cutouts + finder + radar
        cutoutsci = fig.add_subplot(gs[0, 0:5])
        cutouttemp = fig.add_subplot(gs[0, 5:10])
        cutoutdiff = fig.add_subplot(gs[0, 10:15])
        cutoutfinder = fig.add_subplot(gs[0, 15:20])
        radar_ax = fig.add_subplot(gs[0, 25:29], projection="polar")

        # Bottom row: lightcurve + attributes text box
        lc_ax1 = fig.add_subplot(gs[1, 0:20])
        attr_ax = fig.add_subplot(gs[1, 24:32])
        attr_ax.axis("off")

        cutout_fov = _plot_cutouts_and_get_fov(
            cutouts=cutouts,
            axes=[cutoutsci, cutouttemp, cutoutdiff],
        )

        finder_fov_arcsec = _finder_fov_arcsec(cutout_fov, factor=6.0)
        render_finder_stamp(
            cutoutfinder,
            ra,
            dec,
            cache_dir=cutout_cache_dir,
            cache_key=name,
            size=240,
            fov_arcsec=finder_fov_arcsec,
        )

        create_classprob_radar(dict(cls_items), radar_ax, threshold=cls_threshold)

        attr_lines = _make_cls_lines(cls_items, max_items=12)
        attr_lines.append("------------------------")
        attr_lines.extend(info)

        attr_ax.text(
            0.0,
            1.35,
            "\n".join(attr_lines),
            va="top",
            ha="left",
            fontsize=8.5,
            alpha=0.6,
            wrap=True,
        )

    # Layout B: Cutouts + finder stamp only (no classification radar).
    elif has_cutouts and not has_classprobs:
        assert cutouts is not None
        gs = _setup_gs_2x32(fig)

        # Top row: cutouts + finder (wider allocation than Layout A)
        cutoutsci = fig.add_subplot(gs[0, 0:6])
        cutouttemp = fig.add_subplot(gs[0, 6:12])
        cutoutdiff = fig.add_subplot(gs[0, 12:18])
        cutoutfinder = fig.add_subplot(gs[0, 18:24])

        # Bottom row: lightcurve
        lc_ax1 = fig.add_subplot(gs[1, 0:20])

        cutout_fov = _plot_cutouts_and_get_fov(
            cutouts=cutouts,
            axes=[cutoutsci, cutouttemp, cutoutdiff],
        )

        finder_fov_arcsec = _finder_fov_arcsec(cutout_fov, factor=6.0)
        render_finder_stamp(
            cutoutfinder,
            ra,
            dec,
            cache_dir=cutout_cache_dir,
            cache_key=name,
            size=240,
            fov_arcsec=finder_fov_arcsec,
        )

        fig.text(
            0.79,
            0.55,
            "\n".join(info),
            va="top",
            fontsize=8.5,
            alpha=0.6,
        )

    # Layout C: Lightcurve only (no cutouts).
    else:
        lc_ax1 = fig.add_subplot(1, 8, (1, 6))
        fig.subplots_adjust(top=0.8, bottom=0.15)
        plt.subplots_adjust(wspace=0.4, hspace=1.8)

        if has_classprobs:
            info = list(info)
            info.append("------------------------")
            info.append("Class Probabilities:")
            for k, v in cls_items:
                info.append(f"{k}: {v:.3f}")

        fig.text(
            0.795,
            0.75,
            "\n".join(info),
            va="top",
            fontsize=8.5,
            alpha=0.6,
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

    # Give the figure a title (only in simple layout)
    if not has_cutouts:
        fig.suptitle(title if title is not None else f"{name}", fontweight="bold")

    if grid_interval is not None:
        lc_ax1.xaxis.set_major_locator(MultipleLocator(grid_interval))

    lc_ax1.grid(visible=True, axis="both", alpha=0.5)
    lc_ax1.set_ylabel("Magnitude [AB]")
    if not has_cutouts:
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
        tempTab = fluxtable[fluxtable["band"] == fid]
        if len(tempTab) > 0:
            lc_ax1.errorbar(
                tempTab["time"],
                tempTab["mag"],
                tempTab["magerr"],
                color=BANDPASSES[fid]["c"],
                fmt=".",
                markersize=9,
                label=BANDPASSES[fid]["label"],
                mec="black",
                mew=0.5,
            )

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

    # Attribute line (kept close to original behavior)
    if attributes:
        if cutouts:
            fig.text(
                0.5,
                0.985,
                " · ".join(attributes),
                ha="center",
                va="top",
                fontsize=9,
                alpha=0.65,
            )
        else:
            fig.text(
                0.5,
                0.035,
                " - ".join(attributes),
                va="top",
                ha="center",
                fontsize="medium",
                alpha=0.5,
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
    return fig, axes


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
    ny, nx = data2d.shape[-2], data2d.shape[-1]

    try:
        w = WCS(hdr)
        # pixel scales in deg/pix (for each axis)
        scales_deg = proj_plane_pixel_scales(w)  # ndarray, len>=2
        pixscale_arcsec = float(np.mean(scales_deg[:2]) * 3600.0)
        return pixscale_arcsec * float(max(nx, ny))
    except Exception:
        return None


def create_stamp_plot(
    cutouts: Mapping[str, Mapping[str, bytes]], ax, cutout_type: str
) -> float | None:
    data_bytes = next(iter(cutouts.values()))[f"cutout{cutout_type}"]

    data, hdr = load_cutout_image(
        data_bytes,
        return_header=True,
    )

    cutout_fov = cutout_fov_arcsec(data, hdr)

    rotpa = hdr.get("ROTPA", None)
    if rotpa is not None:
        try:
            angle = -float(rotpa)  # close-enough convention: rotate image by -ROTPA
            data = ndimage.rotate(
                data,
                angle=angle,
                reshape=False,  # keep cutout size
                order=1,  # bilinear
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


def render_finder_stamp(
    ax,
    ra: float,
    dec: float,
    *,
    cache_dir: str | None,
    cache_key: str,
    size: int = 240,
    fov_arcsec: float = 45,
) -> None:
    """
    Render a best-effort finder stamp with optional on-disk caching.

    Cache stores a ready-to-imshow 2D array (npz with 'data' + 'label').
    """
    stamp: np.ndarray | None = None
    label: str | None = None
    try:
        if cache_dir is not None:
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
                    stamp2d = normalize_image_for_imshow(stamp)
                    np.savez(cache_path, data=stamp2d, label=label)
                    stamp = stamp2d
        else:
            stamp, label = get_finder_stamp(ra, dec, size=size, fov_arcsec=fov_arcsec)
            if stamp is not None:
                stamp = normalize_image_for_imshow(stamp)
    except Exception:
        stamp, label = None, None

    if stamp is not None:
        ax.imshow(stamp, aspect="equal", cmap="gray")
        add_gap_crosshair(ax, gap_frac=0.2, lw=1.0, alpha=0.9)
        ax.set_title(label or "Finder", fontdict={"fontsize": "small"})
    else:
        ax.text(
            0.5, 0.5, "Finder unavailable", ha="center", va="center", fontsize="small"
        )
        ax.set_title("Finder", fontdict={"fontsize": "small"})

    ax.set_xticks([])
    ax.set_yticks([])


def get_finder_stamp(
    ra: float,
    dec: float,
    size: int = 240,
    surveys: list[str] | None = None,
    timeout: float = 8.0,
    fov_arcsec: float = 45.0,
) -> tuple[np.ndarray | None, str | None]:
    """
    Best-effort finder stamp:
      1) PS1 via get_ps_stamp (PIL image) -> grayscale float array
      2) CDS hips2fits JPEG fallback (DECaLS preferred; SkyMapper/DSS2 later)

    Returns (image_array, label). image_array is 2D float array or None.
    """
    # 1) PS1
    try:
        img = get_ps_stamp(
            ra, dec, size=fov_arcsec / 0.25, color=["y", "g", "i"]
        )  # 0.25 arcsec/px (regarding ztfquery docstring)
        if img is not None:
            arr = np.asarray(img).astype(float)
            if arr.ndim == 3 and arr.shape[2] >= 3:
                arr = 0.2126 * arr[..., 0] + 0.7152 * arr[..., 1] + 0.0722 * arr[..., 2]
            return arr, "PS1"

    except Exception:
        pass

    if surveys is None:
        surveys = [
            "CDS/P/DESI-Legacy-Surveys/DR10/color",
            "CDS/P/DECaLS/DR5/color",
            "CDS/P/DES-DR2/ColorIRG",
            "CDS/P/Skymapper/DR4/color",
            "CDS/P/DSS2/color",
        ]

    labels = {
        "CDS/P/DESI-Legacy-Surveys/DR10/color": "DESI DR10",
        "CDS/P/DES-DR2/ColorIRG": "DES",
        "CDS/P/Skymapper/DR4/color": "SkyMapper",
        "CDS/P/DECaLS/DR5/color": "DECaLS",
        "CDS/P/DSS2/color": "DSS2",
    }

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
            }
            r = requests.get(
                "https://alasky.cds.unistra.fr/hips-image-services/hips2fits",
                params=params,
                timeout=timeout,
            )

            if r.status_code != 200 or not r.content:
                continue

            img = Image.open(io.BytesIO(r.content))
            img = ImageOps.exif_transpose(img).convert("RGB")
            arr = np.asarray(img).astype(float)

            # RGB -> grayscale [0,1]
            arr = 0.2126 * arr[..., 0] + 0.7152 * arr[..., 1] + 0.0722 * arr[..., 2]
            arr /= 255.0

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


def normalize_image_for_imshow(arr: np.ndarray) -> np.ndarray:
    """
    Convert possible image cubes into a 2D (H,W) array suitable for imshow.
    Handles shapes:
      (H,W) -> unchanged
      (H,W,3/4) -> luminance
      (3/4,H,W) -> luminance
      (N,H,W) with N!=3/4 -> take first plane
    """
    a = np.asarray(arr)

    if a.ndim == 2:
        return a

    if a.ndim == 3 and a.shape[2] in (3, 4):
        rgb = a[..., :3].astype(float)
        return 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]

    if a.ndim == 3 and a.shape[0] in (3, 4):
        rgb = a[:3, ...].astype(float)
        return 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]

    if a.ndim == 3:
        return a[0, ...] if a.shape[0] < a.shape[-1] else a[..., 0]

    raise TypeError(f"Unsupported stamp shape for imshow: {a.shape}")


def get_upperlimit_table(dps: Iterable[DataPoint] | None) -> Table | None:
    """
    Build a table of ZTF upper limits from datapoints (id<0 with diffmaglim in body).
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


def clean_classprob_label(label: str) -> str:
    """
    Normalize class probability labels.

    Some classifiers output keys like 'P(SNIa)'. For display purposes we strip
    the outer 'P(...)' if present, leaving 'SNIa'.
    """
    if label.startswith("P(") and label.endswith(")"):
        return label[2:-1]
    return label


def filtered_classprobs(
    classprobs: Mapping[str, float],
    *,
    threshold: float = 0.01,
) -> list[tuple[str, float]]:
    """
    Return (label, prob) pairs with prob >= threshold, sorted descending.
    Labels are cleaned via clean_classprob_label().
    """
    items: list[tuple[str, float]] = []
    for k, v in classprobs.items():
        try:
            fv = float(v)
        except Exception:
            continue
        if fv >= threshold:
            items.append((clean_classprob_label(str(k)), fv))

    items.sort(key=lambda kv: kv[1], reverse=True)
    return items


def create_classprob_radar(
    classprobs: dict[str, float], ax, threshold: float = 0.01
) -> None:
    """
    Radar chart for class probabilities.
    Shows all classes with p >= threshold.
    Gracefully handles N=1/2 via bar plot fallback.
    """
    items = [(k, float(v)) for k, v in classprobs.items() if float(v) >= threshold]
    items.sort(key=lambda kv: kv[1], reverse=True)

    if not items:
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            f"No classes >= {threshold:.2f}",
            ha="center",
            va="center",
            fontsize="small",
        )
        return

    labels = [clean_classprob_label(k) for k, _ in items]

    vals = np.array([v for _, v in items], dtype=float)

    ax.set_ylim(0.0, 1.0)

    rticks = [0.5, 1.0]
    ax.set_yticks(rticks)
    ax.set_yticklabels([str(x) for x in rticks], fontsize=7.5, alpha=0.6)
    ax.tick_params(axis="y", pad=2)

    n = len(vals)
    if n < 3:
        theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
        width = (2 * np.pi) / max(n, 1) * 0.8
        ax.bar(theta, vals, width=width, alpha=0.35, edgecolor="black", linewidth=0.6)
        ax.set_xticks(theta)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_theta_offset(np.pi / 2.0)
        ax.set_theta_direction(-1)
        ax.grid(True, alpha=0.4)
        return

    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    theta_closed = np.r_[theta, theta[0]]
    vals_closed = np.r_[vals, vals[0]]

    ax.set_theta_offset(np.pi / 2.0)
    ax.set_theta_direction(-1)

    ax.plot(theta_closed, vals_closed, linewidth=1.2)
    ax.fill(theta_closed, vals_closed, alpha=0.25)

    ax.set_xticks(theta)
    ax.set_xticklabels(labels, fontsize=8)
    ax.grid(True, alpha=0.4)


class PlotTransientLightcurves(AbsPhotoT3Unit, AbsTabulatedT2Unit):
    """
    Create a (pdf) plot summarizing lightcurves of candidates provided to the unit.

    Features:
      - Optional thumbnails (ZTFCutoutImages / LSSTCutoutImages)
      - Optional finder stamp (PS1 preferred; HiPS fallback)
      - Optional link to TNS (from TNSReports complement)
      - Optional Fritz link (disabled for LSST-only sources)
      - Optional Slack upload
    """

    pdf_path: None | str = None  # Will create random if not set
    titleprefix: str = "AMPEL: "
    save_png: bool = False
    image_cache_dir: str | None = "./image_cache"
    save_dir: str = "./images"
    include_cutouts: bool = False
    fritzlink: bool = True
    webclient: WebClient | None = None
    slack_channel: str | None = None
    slack_token: NamedSecret[str] | None = None

    def post_init(self) -> None:
        os.makedirs(self.save_dir, exist_ok=True)
        if self.image_cache_dir:
            os.makedirs(self.image_cache_dir, exist_ok=True)
        if not self.pdf_path:
            self.pdf_path = tempfile.mkstemp(".pdf", "candidates", self.save_dir)[1]

        if self.slack_channel and self.slack_token is not None:
            self.webclient = WebClient(self.slack_token.get())

    def photz_from_t2s(self, tview: TransientView) -> list:
        """
        Collect information from T2CatalogMatch and T2LSPhotoZ documents, return as list of str.
        """
        photz_list = []
        t2res = tview.get_t2_body(unit="T2DigestRedshifts")
        if isinstance(t2res, dict) and t2res.get("ampel_z", -10) > 0:
            t2cat = tview.get_t2_body(unit="T2CatalogMatch")
            if isinstance(t2cat, dict):
                if (
                    t2cat.get("GLADEv23")
                    and t2cat["GLADEv23"].get("z", None) is not None
                ):
                    photz_list.append("GLADE: {:.2f}".format(t2cat["GLADEv23"]["z"]))
                if t2cat.get("SDSS_spec") and "z" in t2cat["SDSS_spec"]:
                    photz_list.append("SDSS: {:.2f}".format(t2cat["SDSS_spec"]["z"]))
                if t2cat.get("NEDz") and "z" in t2cat["NEDz"]:
                    photz_list.append("NED: {:.2f}".format(t2cat["NEDz"]["z"]))
                if (
                    t2cat.get("wiseScosPhotoz")
                    and "zPhoto_Corr" in t2cat["wiseScosPhotoz"]
                ):
                    photz_list.append(
                        "wiseScos: {:.2f}".format(
                            t2cat["wiseScosPhotoz"]["zPhoto_Corr"]
                        )
                    )
                if t2cat.get("PS1_photoz") and "z_phot" in t2cat["PS1_photoz"]:
                    ps1_photz = float(t2cat["PS1_photoz"]["z_phot"]) / 1000
                    ps1_photzerr = float(t2cat["PS1_photoz"]["z_photErr"]) / 10000
                    photz_list.append(f"PS1: {ps1_photz:.2f} ± {ps1_photzerr:.2f}")
                if t2cat.get("LSPhotoZZou") and "photoz" in t2cat["LSPhotoZZou"]:
                    photz_list.append(
                        "LSZou: {:.2f}".format(t2cat["LSPhotoZZou"]["photoz"])
                    )

            t2cat = tview.get_t2_body(unit="T2LSPhotoZTap")
            if (
                isinstance(t2cat, dict)
                and t2cat.get("T2LSPhotoZTap")
                and "z_phot_median" in t2cat["T2LSPhotoZTap"]
            ):
                photz_list.append(
                    "LS: {:.2f} ± {:.2f}".format(
                        t2cat["T2LSPhotoZTap"]["z_phot_median"],
                        t2cat["T2LSPhotoZTap"]["z_phot_std"],
                    )
                )
        return photz_list

    def attributes_from_t2(
        self,
        tview: TransientView,
        nearby_z: float = 0.02,
        snia_minprob: float = 0.7,
        min_kilonovaness=0,
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
            attributes.append(f"Dist2host: {dist2host:.1f} arcsec")

        t2res = tview.get_t2_body(unit="T2InfantCatalogEval")
        if isinstance(t2res, dict) and t2res.get("action", False):
            attributes.append("InfantEval")

        t2res = tview.get_t2_body(unit="T2RunParsnip")
        if (
            isinstance(t2res, dict)
            and "classification" in t2res
            and t2res["classification"]["SNIa"] > snia_minprob
        ):
            attributes.append("ProbSNIa")

        t2res = tview.get_t2_body(unit="T2KilonovaEval")
        if (
            isinstance(t2res, dict)
            and t2res.get("kilonovaness", -99) > min_kilonovaness
        ):
            attributes.append("Kilonovaness{}".format(t2res["kilonovaness"]))
            attributes.append("LVKmap{}".format(t2res["map_name"]))

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

    def classprobs_from_t2(self, tview: TransientView) -> dict[str, float]:
        """
        Extract class probabilities from T2ClassificationReport (if present).
        """
        t2res = tview.get_t2_body(unit="T2ClassificationReport")
        if not isinstance(t2res, dict):
            return {}

        cls = t2res.get("classification")
        if isinstance(cls, tuple) and len(cls) == 1:
            cls = cls[0]
        if not isinstance(cls, dict):
            return {}

        models = cls.get("models")
        if isinstance(models, (tuple, list)):
            models_iter = models
        elif isinstance(models, dict):
            models_iter = [models]
        else:
            return {}

        out: dict[str, float] = {}
        for m in models_iter:
            if not isinstance(m, dict):
                continue
            probs = m.get("probabilities")
            if not isinstance(probs, dict):
                continue
            for k, v in probs.items():
                try:
                    out[str(k)] = float(v)
                except Exception:
                    continue
        return out

    # Because TransientView does not guarantee any id convention, we introduce these helpers
    # to treat all datapoints as possible photopoints/upperlimits.
    # Eventually, we need to fix that with a proper tabulator-like unit!

    def _iter_t0(self, tview: TransientView) -> Sequence[DataPoint]:
        return tview.t0 or []

    def _is_upperlimit_dp(self, dp: DataPoint) -> bool:
        body = dp.get("body") or {}
        return isinstance(body, dict) and "diffmaglim" in body

    def _get_photopoints_any_id(self, tview: TransientView) -> list[DataPoint]:
        dps = self._iter_t0(tview)
        return [dp for dp in dps if not self._is_upperlimit_dp(dp)]

    def _get_upperlimits_any_id(self, tview: TransientView) -> list[DataPoint]:
        dps = self._iter_t0(tview)
        return [dp for dp in dps if self._is_upperlimit_dp(dp)]

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None
    ) -> UBson | UnitResult:
        buffer = io.BytesIO()
        with PdfPages(self.pdf_path or buffer) as pdf:
            for tran_view in gen:
                photopoints = self._get_photopoints_any_id(tran_view)
                if not photopoints:
                    self.logger.debug("No photopoints", extra={"stock": tran_view.id})
                    continue

                sncosmo_table = self.get_flux_table(photopoints)

                sncosmo_table = sncosmo_table[sncosmo_table["flux"] > 0]
                if sncosmo_table is None or len(sncosmo_table) == 0:
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
                classprobs = self.classprobs_from_t2(tran_view)

                tnsname = None
                if (
                    isinstance(tran_view.extra, dict)
                    and "TNSReports" in tran_view.extra
                ):
                    tnsname = next(iter(tran_view.extra["TNSReports"])).get(
                        "objname", None
                    )

                cutouts = None
                if self.include_cutouts and isinstance(tran_view.extra, dict):
                    for cutout_source in ("LSSTCutoutImages", "ZTFCutoutImages"):
                        raw = tran_view.extra.get(cutout_source)
                        if not isinstance(raw, dict):
                            continue
                        filtered = {k: v for k, v in raw.items() if v}
                        if filtered:
                            cutouts = filtered
                            self.logger.debug(
                                "Found cutouts",
                                extra={"stock": tran_view.id, "source": cutout_source},
                            )
                            break

                fig_from_fluxtable(
                    name,
                    str(ampelid),
                    ra,
                    dec,
                    sncosmo_table,
                    ulim_table,
                    title=title,
                    attributes=attributes,
                    photz_list=photz_list,
                    classprobs=classprobs,
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

        if (
            self.slack_channel is not None
            and self.slack_token is not None
            and self.pdf_path is not None
            and os.path.isfile(self.pdf_path)
        ):
            assert self.webclient is not None

            new_file = self.webclient.files_upload_v2(
                title="My Test Text File", filename="test.pdf", file=self.pdf_path
            )
            fileresp: dict[str, Any] = new_file.get("file", {})
            if len(file_url := fileresp.get("permalink", "")) > 0:
                self.webclient.chat_postMessage(
                    channel=self.slack_channel,
                    text=f"Here is the file: {file_url}",
                )

            try:
                self.webclient.chat_postMessage(
                    channel=self.slack_channel, text="Hello from your app! :tada:"
                )
            except SlackApiError as e:
                assert e.response["error"]

        return None
