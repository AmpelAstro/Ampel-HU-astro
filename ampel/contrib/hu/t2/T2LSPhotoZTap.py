#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2LSPhotoZTap.py
# License           : BSD-3-Clause
# Author            : jnordin
# Date              : 20.04.2021
# Last Modified Date: 21.04.2021
# Last Modified By  : jnordin


from ampel.contrib.hu.t2.T2AbsDatalabMatch import T2AbsDatalabMatch


class T2LSPhotoZTap(T2AbsDatalabMatch):
    """
    Query the NOIR DataLab service for photometric redshifts from the
    Legacy Survey.

    Other queries can in principle be made as long as the string format parameters are the same (ra, dec, match_radius).

    """

    def query(self) -> str:
        return """
            SELECT 
                ra, dec, 
                photo_z.z_phot_median, photo_z.z_phot_mean, photo_z.z_phot_std, 
                photo_z.z_phot_l68, z_phot_u68, photo_z.z_spec, 
                tractor.dered_mag_g, tractor.dered_mag_r, tractor.dered_mag_z, 
                tractor.dered_mag_w1, tractor.dered_mag_w2 , tractor.dered_mag_w3, 
                tractor.dered_mag_w4, tractor.snr_g, tractor.snr_r, tractor.snr_z, 
                tractor.snr_w1, tractor.snr_w2, tractor.snr_w3, tractor.snr_w4 
            FROM 
                ls_dr9.tractor as tractor 
                JOIN ls_dr9.photo_z as photo_z on photo_z.ls_id = tractor.ls_id 
            WHERE 
                't' = Q3C_RADIAL_QUERY(ra, dec,%.6f,%.6f,%.6f)
            """
