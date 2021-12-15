#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2DigestRedshifts.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 06.06.2021
# Last Modified Date: 06.06.2021
# Last Modified By  : jnordin@physik.hu-berlin.de

from typing import List, Dict, Any, Optional, Sequence, Literal, Union
from ampel.struct.UnitResult import UnitResult

from ampel.types import UBson
from ampel.enum.DocumentCode import DocumentCode
#from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.view.T2DocView import T2DocView
from ampel.view.LightCurve import LightCurve
import numpy as np

class T2DigestRedshifts(AbsTiedLightCurveT2Unit):
    """

    Compare potential matches from different T2 units providing redshifts.

    Using table comparisons from (T3) CompareRedshifts to select best matches
    and provide a general uncertainty estimate.

    Available (studied) redshifts are assigned to one of seven redshift "groups", with decreasing average quality.
    The mean redshift of the lowest populated group is returned as 'ampel_z', together with info about this group
    ('group_z_nbr' and 'group_z_precision').

    This is for now an AbsTiedLightCurveT2Unit, but lc info currently not used. Is this needed, or shift to AbsTiedT2Unit?




    """

    # Max redshift uncertainty category: 1-7 (where 7 is any, and 1 only nearby spectroscopic matches)
    max_redshift_category: int


    # Redshift estimates associated with each region ( only rough guideline!!! )
    category_precision: List[float] = [0.0003, 0.003, 0.01, 0.02, 0.04, 0.1, 0.3]


    # CatalogMatch(Local) results might be overriden, for example if specialized catalog is being used
    # Each override dict is assumed to be built asmed to be built according to
    # "catalog_name" : {
    #                    "z_keyword": "redshift field in catalog",
    #                    "max_distance": "max arcsec in which to allow match,
    #                    "max_redshift": "max redshift to use",
    #                    "z_group": "which redshift group to assign to" }
    catalogmatch_override: Optional[Dict[str, Any]]



    # These are the units through which we look for redshifts
    # Which units should this be changed to
    t2_dependency: Sequence[StateT2Dependency[Literal[
		"T2CatalogMatch" ,
		"T2LSPhotoZTap" ,
		"T2CatalogMatchLocal" ,
		"T2MatchBTS"
	]]]




    def _get_lsphotoz_groupz(self, t2_res: Dict[str, any])->List[List[float]]:
        """
        Parse output from T2LSPhotoZTap and investigate whether any matches fulfill group
        redshift criteria.

        Return:
        One list for each of the seven redshift cateogries
        """

        group_z = [[], [], [], [], [], [], []]
        for lsname, lsdata in t2_res.items():
            if lsdata is None:
                continue

            # Warning: all LS checks done with a 10" matching radius, this is thus enforced (in case T2 run with larger radius)
            if lsdata['dist2transient'] > 10:
                self.logger.debug('No Digest redshift LS estimate.', extra={'dist2transient':lsdata['dist2transient']})
                continue

            # First investigate LS spectroscopic redshift
            if lsdata['z_spec'] is not None and lsdata['z_spec'] >- 1:
                if lsdata['z_spec'] < 0.03:
                    # Group I
                    group_z[0].append(lsdata['z_spec'])
                elif lsdata['z_spec'] < 0.15:
                    # Group II
                    group_z[1].append(lsdata['z_spec'])
                elif lsdata['z_spec'] < 0.4:
                    # Group III
                    group_z[2].append(lsdata['z_spec'])
                else:
                    # Group V
                    group_z[4].append(lsdata['z_spec'])
            self.logger.debug('LS debug spec: %s yield %s'%(lsdata, group_z))

            # Now, photometric redshifts
            if lsdata['z_phot_median'] is not None and lsdata['z_phot_median'] >- 1:

                if lsdata['z_phot_median'] < 0.1:
                    # Group IV
                    group_z[3].append(lsdata['z_phot_median'])
                elif lsdata['z_phot_median'] < 0.2:
                    # Group V
                    group_z[4].append(lsdata['z_phot_median'])
                elif lsdata['z_phot_median'] < 0.4:
                    # Group VI
                    group_z[5].append(lsdata['z_phot_median'])
                else:
                    # Group VII
                    group_z[6].append(lsdata['z_phot_median'])
            self.logger.debug('LS debug phot: %s yield %s'%(lsdata, group_z))

        return group_z





    def _get_catalogmatch_groupz(self, t2_res: Dict[str, any])->List[List[float]]:
        """
        Parse output from T2CatalogMatch.

        Made complicated as returns can be both single and lists.


        Return:
        One list for each of the seven redshift cateogries
        """

        group_z = [[], [], [], [], [], [], []]



        for cat_name, cat_matches in t2_res.items():
            if cat_matches is None or cat_matches is False:
                continue
            # Could be list or dict depending on whether the closes or all matches are returned from. To be compatible with both...
            if isinstance(cat_matches, list):
                cat_match_list = cat_matches
            elif isinstance(cat_matches, tuple):
                cat_match_list = list(cat_matches)
            else:
                cat_match_list = [cat_matches]

            for cat_match in cat_match_list:

                # All catalogs have different structure, so doing this individually

                if cat_name == 'NEDz_extcats':
                    if cat_match['dist2transient'] < 2 and cat_match['z'] < 0.03:    # TODO - is this nedz_extcat z ??
                        group_z[0].append(cat_match['z'])
                    elif cat_match['dist2transient'] < 20 and cat_match['z'] < 0.05:
                        group_z[2].append(cat_match['z'])
                    else:
                        group_z[3].append(cat_match['z'])

                if cat_name == 'SDSS_spec':
                    if cat_match['dist2transient'] < 10:    # Implicit restriction as tests where done with this max matching radius
                        group_z[1].append(cat_match['z'])

                if cat_name == 'GLADEv23' and cat_match['dist2transient'] < 10 and cat_match['z'] is not None:    # Implicit restriction as tests where done with this max matching radius
                    if cat_match['z'] < 0.05:
                        group_z[2].append(cat_match['z'])
                    else:
                        group_z[3].append(cat_match['z'])

                if cat_name == 'LSPhotoZZou':
                    # Spec
                    if cat_match['specz'] is not None and cat_match['specz'] >- 0.1:
                        if cat_match['specz'] < 0.15 and cat_match['dist2transient'] < 10:
                            group_z[1].append(cat_match['specz'])
                        elif cat_match['specz'] < 0.2:
                            group_z[2].append(cat_match['specz'])
                        else:
                            group_z[4].append(cat_match['specz'])

                    # Photo-z
                    if cat_match['photoz'] is not None and cat_match['photoz'] >- 0.1:
                        if cat_match['photoz'] < 0.1:
                            group_z[3].append(cat_match['photoz'])
                        elif cat_match['photoz'] < 0.2:
                            group_z[4].append(cat_match['photoz'])
                        elif cat_match['dist2transient'] < 20:
                            group_z[5].append(cat_match['photoz'])
                        else:
                            group_z[6].append(cat_match['photoz'])

                if cat_name == 'wiseScosPhotoz':
                    if cat_match['zPhoto_Corr'] is not None and cat_match['zPhoto_Corr'] >- 0.1:
                        if cat_match['zPhoto_Corr'] < 0.2:
                            group_z[4].append(cat_match['zPhoto_Corr'])
                        else:
                            group_z[5].append(cat_match['zPhoto_Corr'])


                if cat_name == 'twoMPZ':
                    # Photoz
                    if cat_match['zPhoto'] is not None and cat_match['zPhoto'] >- 0.1:
                        if cat_match['zPhoto'] < 0.03:
                            group_z[2].append(cat_match['zPhoto'])
                        else:
                            group_z[3].append(cat_match['zPhoto'])
                    # Specz
                    if cat_match['zSpec'] is not None and cat_match['zSpec'] >- 0.1:
                        group_z[1].append(cat_match['zSpec'])

                if cat_name == 'NEDz' and cat_match['dist2transient'] < 10:    # Implicit restriction as tests where done with this max matching radius
                    if cat_match['z'] < 0.4:
                        group_z[2].append(cat_match['z'])


                # Also check for manual override
                if self.catalogmatch_override:
                    for or_catname, or_catdict in self.catalogmatch_override.items():
                        if or_catname == cat_name:
                            try:
                                cat_z = float(cat_match[or_catdict["z_keyword"]])
                                if float(cat_match['dist2transient']) < or_catdict["max_distance"] and cat_z < or_catdict["max_redshift"]:
                                    group_z[or_catdict["z_group"]-1].append(cat_z)
                            except ValueError:
                                self.logger.info('Cannot parse z', extra={'catdict':cat_match})


        return group_z



    def _get_matchbts_groupz(self, t2_res: Dict[str, any])->List[List[float]]:
        """
        Parse output from T2MatachBTS.

        Any transient with a redshift with two decimals (from SN template matching) is put in Group II,
        those with more (from host, high-res spec) are put into Group I.

        Return:
        One list for each of the seven redshift cateogries
        """

        group_z = [[], [], [], [], [], [], []]

        if 'bts_redshift' in t2_res.keys() and not t2_res['bts_redshift']== '-':
            # BTS redshifts are stored as strings. Crude way to get to redshift precision for evaluation:
            # Take decimal part, remove initial zeroes and cound digits
            decimals = len(t2_res['bts_redshift'].split(".")[1].lstrip("0"))
            z = float(t2_res['bts_redshift'])

            if decimals > 2:
                group_z[0].append(z)
            else:
                group_z[1].append(z)

        self.logger.debug(' bts match yield %s'%(group_z))

        return group_z




    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(self,
        light_curve: LightCurve, t2_views: Sequence[T2DocView]
    ) -> Union[UBson, UnitResult]:
        """

            Parse t2_views from catalogs that were part of the redshift studies.
            Return these together with a "best estimate" - ampel_z

        """

        if not t2_views: # Should not happen actually, T2Processor catches that case
            self.logger.error("Missing tied t2 views")
            return UnitResult(doc_code=DocumentCode.MISSING_INFO)


        # Loop through all potential T2s with redshift information. Each should return an array of arrays, corresponding to redshift maches
        # found in each category. These will be added to sn redshifts
        group_redshifts = [[], [], [], [], [], [], []]


		# Loop through t2_views and collect information.
        for t2_view in t2_views:

            self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
            t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
            # v0.8:
			#t2_res =  t2_view.get_data()

            if t2_view.unit == 'T2LSPhotoZTap':
                new_zs = self._get_lsphotoz_groupz(t2_res)
            elif t2_view.unit == 'T2CatalogMatch' or t2_view.unit == 'T2CatalogMatchLocal':
                new_zs = self._get_catalogmatch_groupz(t2_res)
            elif t2_view.unit == 'T2MatchBTS':
                new_zs = self._get_matchbts_groupz(t2_res)
            else:
                self.logger.error("No instructions for dealing with {}".format(t2_view.unit))
                return UnitResult(doc_code=DocumentCode.MISSING_INFO)

            for k in range(7):
                if len(new_zs[k]) > 0:
                    group_redshifts[k].extend(new_zs[k])
            self.logger.debug('group_z after %s: %s'%(t2_view.unit, group_redshifts))


        # Check for best match
        t2_output = {'group_zs': group_redshifts}
        for k in range(7):
            if (k+1) > self.max_redshift_category:
                # No matches with sufficient precision
                break
            if len(group_redshifts[k]) > 0:
                t2_output['ampel_z'] = np.mean(group_redshifts[k])
                t2_output['group_z_precision'] = self.category_precision[k]
                t2_output['group_z_nbr'] = k+1
                # We then do *not* look for higher group (more uncertain) matches
                break
        if self.catalogmatch_override:
            t2_output["AmpelZ-Warning"]: "Override catalog in use."

        self.logger.debug('digest redshift: %s'%(t2_output))
        return t2_output
