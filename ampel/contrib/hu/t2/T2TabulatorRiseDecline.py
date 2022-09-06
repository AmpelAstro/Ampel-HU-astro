#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2TabulatorRiseDecline.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                28.12.2018
# Last Modified Date:  03.08.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

from typing import Any, Union, TYPE_CHECKING, Sequence, Iterable
import os
import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

from ampel.types import UBson
from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.protocol.LoggerProtocol import LoggerProtocol
from ampel.struct.UnitResult import UnitResult

from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document



def getMag(tab: Table, err=False):
    """
    Shorthand to getting the magnitude from flux.
    """

    m = -2.5* np.log10(tab['flux']) + tab['zp']
    if err:
        # Simple symmetric method
        merr = 2.5 / np.log(10) * (tab['fluxerr']/tab['flux'])
        return m, merr
    else:
        return m

def getMeanflux(tab: Table, jdstart, jdend):
    """
    Return dict with the mean flux for all bands
    with data in the selected range.
    """

    means = {}
    used_bands = []
    tt = tab[ (tab['time'] >= jdstart) & (tab['time'] <= jdend) ]
    for band in set(tt['band']):
        btt = tt[tt['band']==band]
        means[band] = np.average( btt['flux'], weights=1./btt['fluxerr']**2)
        used_bands.append(band)
    return means, used_bands

def getBandBits(bands: Sequence):
    """
    Return number quantifying which bands are included.
    """
    bandval = {'lsstu':1, 'lsstg':2, 'lsstr':4, 'lssti':8, 'lsstz':16, 'lssty':32}
    index = 0
    for band in bands:
        index += bandval[band]
    return index


class T2TabulatorRiseDeclineBase(AmpelBaseModel):
    """
    Derive a number of simple metrics describing the rise, peak and decline of a lc.

    This version assumes input provided by flux table.

    Derived values:
    * t_predetect : time between first detection and previous non-detection
                    in "significant bands". Todo: add min error?
    * t_lc : duration (time between first and most recent detection)
    * jd_max : jd of peak light. None unless bool_peaked
    * jd_det : jd of first detection
    * jd_last : jd of last detection
    * ndet : number of significant detections
    * bool_peaked : is the lc estimated to be declining?
    * bool_pure : has there been no significant non-detections after first detection?
    * bool_rise : was the peak light within "cadence" days of the most recent detection?
   # Not yet implemented    * bool_norise : was the first detection NOT significantly fainter than
        mag_peak IF bool_peaked, ELSE mag_lst
    * bool_hasgaps : The lc has a gap between detections of at least 30 days,
        indicating either a recurrent event or a chance coincidental detection.
    * mag_peak : magnitude at peak light (significant band). Only calculated if bool_peaked
    * mag_det : detection magnitude (significant band)
    * mag_last : magnitude of last detection (significant band)
    Following done for each "color_list" i-ii
    * i-ii_peak : color at peak. None unless bool_peaked AND i+ii obs made within
                 "cadence"  days of jd_max
    * i-ii_det : color at detection. None unless i+ii obs made within "cadence"
                  days of jd_det
    * i-ii_last : color at last detection. None unless i+ii obs made within
                  "cadence" days of jd_last
    * slope_rise_{i,ii} : magnitude slope between jd_det and jd_max. None if bool_norise
    * slope_decline_{i,ii} : magnitude slope between jd_max and jd_lst. None unless bool_peaked

    Additionally, we would like to access host properties like distance and host mags and sizes.
    Would that have to be a different T2?.
    Rather make a different T2 which is chained to the redshift sampler and this.

    "t_cadence" is meant to approximate the rough practical cadence, i.e. how often repeated
    observations should come and what can be seen as "simulateneous"

    "significant_bands" is a list of "deep bands". Collected data from these will be
    used to calculate ages, peak magnitude etc.

    "sigma_det" quantifies what should be seen as a detection.

    "color_list" contains a list of lists of colors which should be considered

    """


    t_cadence: float = 5.
    significant_bands: Sequence[str] = ['lsstg', 'lsstr', 'lssti', 'lsstz']
    sigma_det: float = 5.
    sigma_slope: float = 3. # Threshold for having detected a slope
    color_list: Sequence[Sequence[str]] = [['lsstu','lsstg'], ['lsstg','lsstr'],
                        ['lsstr','lssti'], ['lssti','lsstz'],['lsstz','lssty']]
    max_tgap: int = 30


    if TYPE_CHECKING:
        logger: LoggerProtocol

    def get_bandfeatures(self, ftable):
        """
        Go through all bands of input table and:
        - Try to find a peak, together with rise and fall slopes.
        - Based on this try to determine if it is rising, falling and have peaked.
        - Storing slopes where possible.
        """

        def linearFunc(x,intercept,slope):
            y = intercept + slope * x
            return y

        banddata = {}
        tscale = np.mean(ftable['time'])

        for band in set(ftable['band']):
            bt = ftable[ftable['band'] == band]

            max_flux = bt['flux'].max()
            max_flux_time = bt[bt['flux']==max_flux]['time'][0]
            banddata['jd_peak_'+band] = max_flux_time

            # Divide Table
            riset = bt[ bt['time']<=max_flux_time ]
            fallt = bt[ bt['time']>=max_flux_time ]

            # Examine rise
            if len(riset)>1:
                fit, cov=curve_fit(linearFunc,riset['time']-tscale,riset['flux'],
                                    sigma=riset['fluxerr'],absolute_sigma=True)
                banddata['rise_slope_'+band] = fit[1]
                banddata['rise_slopesig_'+band] = fit[1] / np.sqrt(cov[1][1])
            if len(fallt)>1:
                fit, cov=curve_fit(linearFunc,fallt['time']-tscale,fallt['flux'],
                                    sigma=fallt['fluxerr'],absolute_sigma=True)
                banddata['fall_slope_'+band] = fit[1]
                banddata['fall_slopesig_'+band] = fit[1] / np.sqrt(cov[1][1])

        # In v1 we also had a requirement that a sufficient time should have
        # passed after peak, such that we "should" have seen a decline.
        # We here skip this, hopefully the slope fit is sufficient

        # Check whether we have a significant rise detected in any band.
        risepulls = [ banddata.get('rise_slopesig_'+band,0)
                        for band in set(ftable['band']) ]
        if sum(risepulls)>self.sigma_slope:
            banddata['bool_rise'] = True
        else:
            banddata['bool_rise'] = False
        # Check whether we see a decline
        decpulls = [ banddata.get('fall_slopesig_'+band,0)
                        for band in set(ftable['band']) ]
        if sum(decpulls) < -self.sigma_slope:
            banddata['bool_fall'] = True
        else:
            banddata['bool_fall'] = False

        # If the transient has both a rise and a fall we can
        # define a central peak
        if banddata['bool_rise'] and banddata['bool_fall']:
            banddata['bool_peaked'] = True
            # Use jd of all bands for which we could estimate rise+fall
            banddata['jd_peak'] = np.median(
                                [ banddata['jd_peak_'+band]
                                       for band in set(ftable['band'])
                                  if 'rise_slope_'+band in banddata and
                                     'fall_slope_'+band in banddata] )
        else:
            banddata['bool_peaked'] = False

        # Could include the abs mag at peak, but we argued this would not
        # be as useful?

        return banddata



    def compute_stats(self, flux_table: Table) -> dict[str, Any]:

        # Output dict that we will start to populate
        o: dict[str, Any] = {}

        # Step 1. Base determinations based on combined detections
        self.logger.debug("Starting joint band RiseDeclineStat estimations")

        # Create subset of table with significant detections in significant bands
        band_mask = [bandobs in self.significant_bands for bandobs in flux_table['band'] ]
        sig_mask = (flux_table['flux']) / flux_table['fluxerr']>self.sigma_det
        det_table = flux_table[band_mask & sig_mask]

        o['ndet'] = len(det_table)
        if o['ndet']==0:
            o['success'] = False
            o['cause'] = "No data survive significance criteria."
            return o

        o["jd_det"] = det_table['time'].min()
        o["jd_last"] = det_table['time'].max()

        # Get the max time of obs in signifant bands prior to jd_det
        if flux_table[band_mask]['time'].min()<o['jd_det']:
            o["t_predetect"] = o['jd_det'] - flux_table[band_mask][
                    flux_table['time'][band_mask]<o["jd_det"] ]['time'].max()
        else:
            o["t_predetect"] = None

        o["mag_det"] = float( getMag(flux_table[flux_table['time']==o["jd_det"]]) )
        o["band_det_id"] = getBandBits( [flux_table[flux_table['time']==o["jd_det"]]['band'][0]] )

        o["mag_last"] = float( getMag(flux_table[flux_table['time']==o["jd_last"]]) )
        o["band_last_id"] = getBandBits( [flux_table[flux_table['time']==o["jd_last"]]['band'][0]] )

        o["t_lc"] = o["jd_last"] - o["jd_det"]

        # Check for non-signficant obs between det and last
        ultab = flux_table[band_mask & ~sig_mask]
        if sum( (ultab['time']>=o["jd_det"]) & (ultab['time']<=o["jd_last"]) )==0:
            # No non-detection among signifcant bands between start and end
            o["bool_pure"] = True
        else:
            o["bool_pure"] = False

        # We start measuring bandfeatures t_cadence days prior to first detection
        # to allow some nondetection data to be included.
        # (Note: we previously only used significant bands here. Prob wrong. )
        time_mask = ( (flux_table['time']>(o['jd_det']-self.t_cadence)) &
                      (flux_table['time']<(o['jd_last']+self.t_cadence)) )
        o.update( self.get_bandfeatures(flux_table[time_mask]) )


        # If there is a peak we additionally check whether this is within t_cadence
        # days of detector or last, and call this fastrise and fastfall
        o['bool_fastrise'], o['bool_fastfall'] = None, None
        if o['bool_peaked']:
            if np.abs(o['jd_det']-o['jd_peak'])<self.t_cadence:
                o['bool_fastrise'] = True
            else:
                o['bool_fastrise'] = False
            if np.abs(o['jd_last']-o['jd_peak'])<self.t_cadence:
                o['bool_fastfall'] = True
            else:
                o['bool_fastfall'] = False
            o['t_rise'] = o['jd_peak'] - o['jd_det']
            o['t_fall'] = o['jd_last'] - o['jd_peak']

        # Are there long gaps among the detections?
        jdsorted = np.unique(flux_table["time"])
        if len(jdsorted) > 1:
            if (jdsorted[1:] - jdsorted[0:-1]).max() > self.max_tgap:
                o["bool_hasgaps"] = True
            else:
                o["bool_hasgaps"] = False
        else:
            o["bool_hasgaps"] = None

        # Color
        # Define time subsets at detection, last (significant) and peak (if defined)
        # In each, get the mean flux in each band.
        # We assume that the zeropoint is the same for these fluxes!
        # Also, only exist for positive fluxes (...)

        # Detection colors
        fluxdict, fluxbands =  getMeanflux(flux_table, o["jd_det"]-self.t_cadence/2,
                                            o["jd_det"]+self.t_cadence/2)
        o['det_bands'] = getBandBits(fluxbands)
        for colbands in self.color_list:
            if fluxdict.get(colbands[0],-1)>0 and fluxdict.get(colbands[1],-1)>0:
                o[f"{colbands[0]}-{colbands[1]}_det"] = -2.5 * np.log10(fluxdict[colbands[0]] / fluxdict[colbands[1]])
        # Last obs colors
        fluxdict, fluxbands =  getMeanflux(flux_table, o["jd_last"]-self.t_cadence/2,
                                            o["jd_last"]+self.t_cadence/2)
        o['last_bands'] = getBandBits(fluxbands)
        for colbands in self.color_list:
            if fluxdict.get(colbands[0],-1)>0 and fluxdict.get(colbands[1],-1)>0:
                o[f"{colbands[0]}-{colbands[1]}_last"] = -2.5 * np.log10(fluxdict[colbands[0]] / fluxdict[colbands[1]])
        # Peak colors, if found
        if o['bool_peaked']:
            fluxdict, fluxbands =  getMeanflux(flux_table, o["jd_peak"]-self.t_cadence/2,
                                                o["jd_peak"]+self.t_cadence/2)
            o['peak_bands'] = getBandBits(fluxbands)
            for colbands in self.color_list:
                if fluxdict.get(colbands[0],-1)>0 and fluxdict.get(colbands[1],-1)>0:
                    o[f"{colbands[0]}-{colbands[1]}_peak"] = -2.5 * np.log10(fluxdict[colbands[0]] / fluxdict[colbands[1]])


        o["success"] = True
        return o


class T2TabulatorRiseDecline(AbsStateT2Unit, AbsTabulatedT2Unit, T2TabulatorRiseDeclineBase):

    plot_prob: float = 0.
    path_testplot: str = "/home/jnordin/tmp/t2test/"

    def test_plot(self, name, table, t2result):
        """
        for debugging

        Create panel for each band, showing data + the significant dets.

        Mark times for detection, last det + peak (if set)
        Somehow indicate tilt (from peak pos?)
        Write summary of boolean conclusion.

        Save as stockid + ndet

        """

        import matplotlib.pyplot as plt

        bands = set(table['band'])

        fig, axs = plt.subplots(2, 3)

        for k, band in enumerate(bands):
            bt = table[table['band']==band]
            ax = axs[int(k/3), k%3]

            ax.errorbar(bt['time'], bt['flux'], yerr=bt['fluxerr'],
                            fmt='o', label=band)

            # Detection times
            if 'jd_det' in t2result:
                ax.axvline(t2result["jd_det"], color='grey', linewidth=3, alpha=0.5 )
            if 'jd_last' in t2result:
                ax.axvline(t2result["jd_last"], color='grey', linewidth=3, alpha=0.5 )
            if t2result['bool_peaked']:
                ax.axvline(t2result["jd_peak"], color='red', linewidth=2, alpha=0.8 )
            if 'jd_peak_'+band in t2result:
                ax.axvline(t2result["jd_peak_"+band], color='green', linewidth=1)

            # We next wich to indicate the slopes we have measured
            if 'rise_slope_'+band in t2result:
                slope = t2result['rise_slope_'+band]
                startpoint = t2result['jd_det']
                xvals = np.array([-t2result['t_lc']/4, t2result['t_lc']/4])
                yvals = xvals * t2result['rise_slope_'+band] + np.mean(bt['flux'])
                xvals += (t2result['jd_det']+t2result['t_lc']/4)
                col = 'black'
                if np.abs(t2result['rise_slopesig_'+band])>3:
                    col = 'red'
                ax.plot(xvals,yvals,color=col)
            if 'fall_slope_'+band in t2result:
                xvals = np.array([-t2result['t_lc']/4, t2result['t_lc']/4])
                yvals = xvals * t2result['fall_slope_'+band] + np.mean(bt['flux'])
                xvals += (t2result['jd_last']-t2result['t_lc']/4)
                col = 'black'
                if np.abs(t2result['fall_slopesig_'+band])>3:
                    col = 'red'
                ax.plot(xvals,yvals,color=col)


#            ax.legend()
            ax.set_xlabel('MJD')
            ax.set_ylabel(band)

            # Create text string
            title = "ndet: %s " % (
                t2result["ndet"]
                )

            for boolprop in ["peaked", "pure", "rise", "hasgaps", "fall", "fastrise", "fastfall"]:
                if t2result[f"bool_{boolprop}"]:
                    title += f"{boolprop} "

            ax.set_title(title,{'fontsize':8})

        # Store figure
        path = os.path.join(self.path_testplot, "{}_{}.pdf".format(name, t2result['ndet']))
        plt.tight_layout()
        plt.savefig(path)
        plt.clf()



    def process(self,
        compound: T1Document,
        datapoints: Iterable[DataPoint],
        ) -> Union[UBson, UnitResult]:
        """
        Process datapoints belonging to one state of one transient.
        A commong Table is generated which is used as input
        to the feature generator.
        """

        # Convert input datapoints to standardized Astropy Table
        # Using standard tabulators
        flux_table = self.get_flux_table(datapoints)


        # Calculate get_features
        features = self.compute_stats(flux_table)

        #if self.do_testplot:
        if features['success'] and np.random.uniform()<self.plot_prob:
            self.test_plot(compound.get('stock'), flux_table, features)

        return features
