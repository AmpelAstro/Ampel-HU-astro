#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2GetLensSNParameters.py
# License           : BSD-3-Clause
# Author            : alice.townsend@physik.hu-berlin.de
# Date              : 22.11.2021
# Last Modified Date: 19.10.2022
# Last Modified By  : alice.townsend@physik.hu-berlin.de


import numpy as np
from astropy.table import Table
from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo
from bisect import bisect_left

class T2GetLensSNParameters(T2RunSncosmo):
    unc: float

    def _get_fit_metrics(self, sncosmo_result, sncosmo_table, sncosmo_model) -> dict:
        #-------------------------------- Overriden method -----------------------------#
        lc_metrics = super()._get_fit_metrics(sncosmo_result, sncosmo_table, sncosmo_model)
        # Record the closest detection to peak in each band
        t0 = sncosmo_model.get('t0') # get time at peak of light curve
        sncosmo_table.add_column('temporary', name='epoch', index = 0) # add column with temp values for epoch
        sncosmo_table.sort(['time']) # sort data chronologically
        
        def take_closest(myList, myNumber):
            '''
            Assumes myList is sorted. Returns closest value to myNumber.
            If two numbers are equally close, return the smallest number.
            '''
            pos = bisect_left(myList, myNumber)
            if pos == 0:
                return myList[0]
            if pos == len(myList):
                return myList[-1]
            before = myList[pos - 1]
            after = myList[pos]
            if after - myNumber < myNumber - before:
                return after
            else:
                return before
        
        for band in np.unique(sncosmo_table['band']):
            table_name = str(band) + 'cut'
            globals()[table_name] = Table(sncosmo_table[0:0]) # create empty table for each band
            band_table = sncosmo_table[(sncosmo_table['band'] == band)] # split data into different bands
            time_minus7 = take_closest(band_table['time'], t0 - 7) # find closest time to t0-7 days in each band
            if (t0 - 7 - self.unc) <= time_minus7 <= (t0 - 7 + self.unc):
                band_table['epoch'][np.where(band_table['time'] == time_minus7)] = 'minus7'
                for row in band_table[np.where(band_table['time'] == time_minus7)]:  # iterate over temporary table that is a subset of data  
                    globals()[table_name].add_row(row)
            time_closest = take_closest(band_table['time'], t0) # find closest time to peak in each band
            if (t0 - self.unc) <= time_closest <= (t0 + self.unc):
                band_table['epoch'][np.where(band_table['time'] == time_closest)] = 't0'
                for row in band_table[np.where(band_table['time'] == time_closest)]:  # iterate over temporary table that is a subset of data  
                    globals()[table_name].add_row(row)
            time_plus7 = take_closest(band_table['time'], t0 + 7) # find closest time to t0+7 days in each band
            if (t0 + 7 - self.unc) <= time_plus7 <= (t0 + 7 + self.unc):
                band_table['epoch'][np.where(band_table['time'] == time_plus7)] = 'plus7'
                for row in band_table[np.where(band_table['time'] == time_plus7)]:  # iterate over temporary table that is a subset of data  
                    globals()[table_name].add_row(row)

        # Calculate the colour at different epochs
        def calculate_colour_peak(band1, band2, epoch):
            '''
            Calculate the colour from two different bands of a user's choosing.
            Bluer colour is band1.
            epoch is the epoch within which we are searching (i.e. t0, minus7, or plus7).
            '''
            band1_table = str(band1) + 'cut'
            band2_table = str(band2) + 'cut'
            if band1 not in sncosmo_table['band'] or band2 not in sncosmo_table['band']:
                return None, None
            elif epoch not in globals()[band1_table]['epoch'] or epoch not in globals()[band2_table]['epoch']:
                return None, None
            else:
                colour = -2.5 * np.log10(globals()[band1_table][np.where(globals()[band1_table]['epoch'] == epoch)]['flux']
                                 / globals()[band2_table][np.where(globals()[band2_table]['epoch'] == epoch)]['flux'])
                unc1 = 2.5 * 0.434 * (globals()[band1_table][np.where(globals()[band1_table]['epoch'] == epoch)]['fluxerr']
                              /globals()[band1_table][np.where(globals()[band1_table]['epoch'] == epoch)]['flux'])
                unc2 = 2.5 * 0.434 * (globals()[band2_table][np.where(globals()[band2_table]['epoch'] == epoch)]['fluxerr']
                              /globals()[band2_table][np.where(globals()[band2_table]['epoch'] == epoch)]['flux'])
                unc = (unc1**2 + unc2**2)**0.5
                return colour.data[0], unc.data[0]
        
        lc_metrics['r_i_colour_peak'], lc_metrics['r_i_colour_peak_err'] = calculate_colour_peak('ztfr', 'ztfi', 't0')
        lc_metrics['r_i_colour_plus7'], lc_metrics['r_i_colour_plus7_err'] = calculate_colour_peak('ztfr', 'ztfi', 'plus7')
        lc_metrics['r_i_colour_minus7'], lc_metrics['r_i_colour_minus7_err'] = calculate_colour_peak('ztfr', 'ztfi', 'minus7')
        lc_metrics['g_r_colour_peak'], lc_metrics['g_r_colour_peak_err']  = calculate_colour_peak('ztfg', 'ztfr', 't0')
        lc_metrics['g_r_colour_plus7'], lc_metrics['g_r_colour_plus7_err'] = calculate_colour_peak('ztfg', 'ztfr', 'plus7')
        lc_metrics['g_r_colour_minus7'], lc_metrics['g_r_colour_minus7_err'] = calculate_colour_peak('ztfg', 'ztfr', 'minus7') 
        lc_metrics['g_i_colour_peak'], lc_metrics['g_i_colour_peak_err']  = calculate_colour_peak('ztfg', 'ztfi', 't0')
        lc_metrics['g_i_colour_plus7'], lc_metrics['g_i_colour_plus7_err'] = calculate_colour_peak('ztfg', 'ztfi', 'plus7')
        lc_metrics['g_i_colour_minus7'], lc_metrics['g_i_colour_minus7_err'] = calculate_colour_peak('ztfg', 'ztfi', 'minus7')        

        #Calculate observed magnitude close to peak
        def calculate_obsmag_peak(band1, epoch):
            '''
            Calculates observed magnitude for a particular epoch
            '''
            band1_table = str(band1) + 'cut'
            if band1 not in sncosmo_table['band']:
                return None
            elif epoch not in globals()[band1_table]['epoch']:
                return None
            else:
                obs_mag = -2.5 *np.log10( globals()[band1_table][np.where(globals()[band1_table]['epoch'] == epoch)]['flux'] ) + 25
                return obs_mag.data[0]
        
        lc_metrics['obsmag_ztfg_peak'] = calculate_obsmag_peak('ztfg', 't0')
        lc_metrics['obsmag_ztfr_peak'] = calculate_obsmag_peak('ztfr', 't0')
        lc_metrics['obsmag_ztfi_peak'] = calculate_obsmag_peak('ztfi', 't0')
        lc_metrics['obsmag_ztfg_plus7'] = calculate_obsmag_peak('ztfg', 'plus7')
        lc_metrics['obsmag_ztfr_plus7'] = calculate_obsmag_peak('ztfr', 'plus7')
        lc_metrics['obsmag_ztfi_plus7'] = calculate_obsmag_peak('ztfi', 'plus7')
        lc_metrics['obsmag_ztfg_minus7'] = calculate_obsmag_peak('ztfg', 'minus7')
        lc_metrics['obsmag_ztfr_minus7'] = calculate_obsmag_peak('ztfr', 'minus7')
        lc_metrics['obsmag_ztfi_minus7'] = calculate_obsmag_peak('ztfi', 'minus7')

        return lc_metrics