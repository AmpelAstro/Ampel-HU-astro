#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/DevSkyportalWebhookPublisher.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                17.02.2023
# Last Modified Date:  14.03.2023
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from itertools import islice
from typing import Iterable, TYPE_CHECKING, Literal, Any, Optional
from collections.abc import Generator
import datetime
from astropy.time import Time
import json
import base64
import requests
import joblib, io

from ampel.types import ChannelId, UBson, T3Send
from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.struct.UnitResult import UnitResult
from ampel.struct.StockAttributes import StockAttributes
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.enum.DocumentCode import DocumentCode
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.secret.NamedSecret import NamedSecret
from ampel.view.TransientView import TransientView
from ampel.view.T2DocView import T2DocView
from ampel.util.mappings import get_by_path
from ampel.util.compression import decompress
from ampel.struct.T3Store import T3Store
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper

if TYPE_CHECKING:
    from ampel.content.JournalRecord import JournalRecord


class SelectT2Skyportal(AmpelBaseModel):
    """
    Define which T2 documents and fields therein to propagate.
    """
    publish_plots: bool = False
    config: None | int = None
    jd_end: None | float = None
    results_mapping: dict[str,Any] = {}
    message_keys: list[str]

class DevSkyportalWebhookPublisher(AbsPhotoT3Unit):
    """

    Dev unit for assembling information for publishing to the Skyportal
    webhook endpoint.
    
    We assume that there is a direct mapping between a t2unit config running
    
    See:
    https://skyportal.io/docs/analysis.html#external-analysis-services
    https://github.com/skyportal/skyportal/pull/3918
    
    Expected format (Feb 17)
    params={
  	"show_parameters": True,     # What does this do?
  	"analysis": {
    		"plots": ...,
    		"inference_data": ...,
    		"results": {
      			"format": "json",
      			"data": {
        			"external_provenance": {
          				"hash": "23baef56",
          				"spectrum_ids": [34, 68, 125],
          				"model": "jsbFitter23"
        			},
        		"classification": "DY Per",
        		"classification_proba": 0.94
      			}
    		}
  	}
    }

    url = f"http://<url>:5000/api/obj/MySourceName/analysis_upload/{service_id}"
    r = requests.post(url, headers={'Authorization':   'token uuid'}, json=params)

    * Our responsibility to ensure provenance, deletion of old analysis and rights propagation? *

      
    AMPEL provenance:
    external_provenance should include:
    - origin (AMPEL)
    - Name of unit
    - Config hash
    - link
    - stock
    - channel
    - code 
    - last meta timestamp
    
    Selecting what to upload
    t2_publish = { 'unitname': {
    	publish_plots: bool = False, # push any entries stored under "plots",
    	config: int | None = None,    # the 
    	results_mapping: {
    	    'output_key': ['path','to','val'],
        }
    }

    Q:
    - How to determine the object id? Currently assumes ZTFIdMapper maps to the skyportal source name.
    - Should we require a specific config to be given as parameter, or is the {channel,unit} combo sufficient?
          i.e. Is it possible that a channel might run multiple copies of a T2Unit, but only want the result of 
               one to be pushed? I guess that could be (think different redshift priors etc). 
               
    Answered Q:
    - How to keep track of the service_id? (fixed list maintained in this file? configuration?)
      Cannot keep them in the file, anyone can generate thos. But also does not need to be a Secret.
      So keep as standard unit parameter.
    - Could the output of more than one T2 be included as on analysis service? 
      Probably too complex, would need to build a new config etc. So we will have one Skyportal call per t2.
    - How should be verify whether a new analysis should be uploaded, or possible the old removed? 
      For now we will simply delete any existing analysis with the same same analysis service id, as only 
      very few are allowed concurrently. Could change 
      

    """

    # Identification
    access_token: str = 'token'    # Will be a secret
    service_id: int        
    skyportal_url: str = 'url'
    provenance_origin: str = 'AMPEL' # Will be propagated as part of provenance. Possibly replace with e.g. 'AMPEL@DESY' 

    # Need some mechanism for determining MySourceName
    # name_filter: dict[str, str] = {'ZTF name': 'ZTF', 'TNS ID': 'TNS'}
 
    # Mapping of t2units and contents to parse for submission data. 
    t2_publish: dict[str, SelectT2Skyportal]

    # Dev
    debug: bool = False         # Only print/store, no publish    
    
    
    def get_external_provenance(self, t2view: T2DocView) -> dict:
        """
        For a transient view, construct a dict which contains 
        the necessary information for skyportal analysis service provenance:
            {
            	'origin': 'AMPEL',
            	'unit': 'T2RunParsnip',
            	'config': hash,
            	'link': link,
            	'stock': stock,
            	't2_type': t2_type         
            	('channel': channel,     # Is this needed? Not included in this run)
            	'code': 0, # Maybe obvious?
            	'meta_ts':  timestamp # latest timestamp in meta
            }
        Could this release more information than we wish? I guess this comes back to only pushing 
        to groups with the access you want.        
        
        TODO: We assume config is an integer hash id, which we turn into a string.
        """
        prov = {'origin': self.provenance_origin, 'unit': t2view.unit,
                'stock': t2view.stock, 'link':t2view.link, 'code': t2view.code, 't2_type': t2view.t2_type}
        # More complex - for config we would like the hash rather then full dict
        print('TODO rehash config if needed?')
        prov['config'] = t2view.config
        
        prov['ts'] = t2view.get_time_updated()
        
        return prov
        
    def _get_skyportal_analysis(self, source_name: str) -> list | None:
        """
          Search for existing analysis results for this source and the service id associated 
          to this unit.
          Will _not_ match t2 config, but assume that service id is associated only to one.
        """
        url = f"http://{self.skyportal_url}/api/sources/{source_name}"
        r = requests.get(url, headers={'Authorization':   'token '+self.access_token}, params={"includeAnalyses" : True})
        if not r.ok:
            # Either no connection or source not existing
            return None
        skyportal_analysis = [a for a in r.json()['data']['analyses'] if int(a['analysis_service_id'])==self.service_id]
        return skyportal_analysis

    def _del_skyportal_analysis(self, skyportal_analysis: list) -> None:
        """
            Delete existing skyportal analysis part of the list.
        """
        for an in skyportal_analysis:
            url = f"http://{self.skyportal_url}/api/obj/analysis/{an['id']}"
            print('deleting this', url)
            r = requests.delete(url, headers={'Authorization':   'token '+self.access_token})
            
        return None


    def _post_skyportal_analysis(self, source_name: str, analysis_payload: dict) -> int | None:
        """
            Publish the payload as an analysis to a Skyportal instance 
        """
            
        # Post payload, will go something like (but through client class?)
        service_url = f"http://{self.skyportal_url}/api/obj/{source_name}/analysis_upload/{self.service_id}"
        r = requests.post(service_url, headers={'Authorization':   'token '+self.access_token}, json=analysis_payload)
        if r.ok:
            return r.json()['data']['id']
        else:
            return None

        
        
        
    def process(self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None) -> UBson | UnitResult:
        """
        Loop through provided TransientViews, extract data and publish according to the
        configured schema.
        """
        
        for tran_view in gen:

            # Chategorized by stock (internal ampel ID) and channel
            assert tran_view.stock is not None
            stock = tran_view.id
            source_name = ZTFIdMapper.to_ext_id(stock)
            print('Assumed skyportal source name', source_name)
            # TODO: Have the mapper as unit variable?
            
            # Get analysis id and submission times from a skyportal instance
            skyportal_info = self._get_skyportal_analysis(source_name)
            # TODO: If a source does not exist in the skyportal, publish this first? 
            if skyportal_info is None:
                print('Source not found - try to post this first?')
                print('send error message to T3Send')
                continue
                        
            # Should channels be published?
            channels = tran_view.stock.get("channel")
            assert channels is not None
            channel = channels[0] if isinstance(channels, (list, tuple)) and len(channels) == 1 else '/'.join([str(c) for c in channels])

            
            for t2unit, publish_settings in self.t2_publish.items():
                payload: dict[str, Any] = {'show_parameters': True, 'analysis': {'results':{'format':'json','data':{}}}}
                for t2view in tran_view.get_t2_views(unit=t2unit):

                    # So, what if we have mutiple copies of the same t2unit run for this channel, but only one config
                    # should be uploaded
                    print('TODO: Check whether this is the correct t2 version.')
                               
                    # Get the provenance info
                    ext_prov = self.get_external_provenance(t2view)
                    print('DEBUG for {} got provenance {}'.format(t2unit, ext_prov))
                    if not ext_prov['code']==0:
                        print('... t2 not completed')
                        continue
                        
                    # Require config matching if set
                    if publish_settings.config and not publish_settings.config==ext_prov['config']:
                        print('DEBUG - config not matching, skipping.')
                        continue
                        
                    # If ext_prov is already set, something has gone wrong (should only get here once)
                    if 'external_provenance' in payload['analysis']['results']['data']:
                        # We assume that we are here for correct reasons (wanted configuration)
                        # and we will here only check for the latest timestamp
                        if ext_prov['ts']>payload['analysis']['results']['data']['external_provenance']['ts']:
                            print('... found later result, using this')
                            payload['analysis']['results']['data']['external_provenance'] = ext_prov
                        else:
                            print('... this doc made earlier, ignoring.')
                            continue
                    else:
                        payload['analysis']['results']['data']['external_provenance'] = ext_prov


                    # The below checks could become relevant if we do not automatically delete any old analysis with 
                    # the same analysis id. As this is the case now, these checks can prob be skipped.
                                    
                    # Have some results, check what already is uploaded 
                    print('TODO: Check skyportal for this {object,service}, get our timestap of latest submission')
                
                    # Check whether the latest timestamp of this T2Record is later than what has been published 
                    print('TODO: Published timestamp older than we have, otherwise we would have continued here.')
                    db_ts = t2view.get_time_updated()
                    print('t2 timestamp', db_ts)
                    print('existing skyportal modified analysis timestamps', [
                        int(datetime.datetime.fromisoformat(a['modified']).timestamp()) for a in skyportal_info] )
                    
                    # Collect other info - plots and results
                    body = t2view.get_payload()
                    if body is None:
                        continue 
                        
                    # Plots
                    if publish_settings.publish_plots and (plotlist:=body.get('plots')):
                        payload['analysis']['plots'] = []
                        for plot in plotlist:
                            raw_im = decompress(plot['svg'])
                            payload['analysis']['plots'].append({"format": "svg+xml", "data": base64.b64encode(raw_im).decode()})

		    # Add analysis information. For now we combine all of them also to the message.                     
                    message = ''
                    for label, path in publish_settings.results_mapping.items():
                        if result := get_by_path(body, path):
                            payload['analysis']['results']['data'][label] = result
                            if label in publish_settings.message_keys:
                                if type(result)==float:
                                     message = message+f'{label}: {float(f"{result:.2g}"):g}. '
                                elif type(result)==int or type(result)==str:
                                    message = message+'{}: {} '.format(label,result)
                    payload["message"] = message

                # Verify that something was found
                if not 'external_provenance' in payload['analysis']['results']['data']:
                    print('DEBUG - nothing found for this unit.')
                    continue

                            
                # Payload for this unit ready.
                if self.debug:
                    print('Payload:')
                    payload['analysis']['plots'][0]['data'] = payload['analysis']['plots'][0]['data'][0:20]
                    print(payload)
                else:
                    # Ready to upload a new analysis. First remove any existing:
                    self._del_skyportal_analysis(skyportal_info)
                    # Upload new
                    analysis_id = self._post_skyportal_analysis(source_name, payload)
                    print('got new analysis id, should be _sent_ to the object journal', analysis_id)
                     


        return None
