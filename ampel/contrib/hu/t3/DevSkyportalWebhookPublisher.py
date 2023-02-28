#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/DevSkyportalWebhookPublisher.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                17.02.2023
# Last Modified Date:  17.02.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from itertools import islice
from typing import Iterable, TYPE_CHECKING, Literal, Any, Optional
from collections.abc import Generator
import datetime
from astropy.time import Time
import json
import base64


from ampel.types import ChannelId, UBson, T3Send
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

if TYPE_CHECKING:
    from ampel.content.JournalRecord import JournalRecord


class DevSkyportalWebhookPublisher(AbsPhotoT3Unit):
    """

    Dev unit for assembling information for publishing to the Skyportal
    webhook endpoint.
    
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

    Notes:
    - Get an analysis_id (same as service_id) through api/analysis_service
    - analysis expect to contain
      * interference_data (`arviz`). Looks like we can leave this out.
      * plots: list of plots
      * "freeform python object" packaged with joblib.dump
      * all three of these data are encoded with base64.b64encode (are they already at our end?)
      
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
    	publish_plots: bool, # push any entries stored under "plots"
    	results_mapping: {
    	    'output_key': ['path','to','val'],
        }
    }

    Q:
    - How should be verify whether a new analysis should be uploaded, or possible the old removed?
    - How to determine the object id? (presumably some name mapper)
    - Should we require a specific config to be given as parameter, or is the {channel,unit} combo sufficient?
          i.e. Is it possible that a channel might run multiple copies of a T2Unit, but only want the result of 
               one to be pushed? I guess that could be (think different redshift priors etc). 
               
    Answered:
    - How to keep track of the service_id? (fixed list maintained in this file? configuration?)
      Cannot keep them in the file, anyone can generate thos. But also does not need to be a Secret.
      So keep as standard unit parameter.
    - Could the output of more than one T2 be included as on analysis service? 
      Probably too complex, would need to build a new config etc. So we will have one Skyportal call per t2.

    """

    # Identification
    access_token: str = 'token'    # Will be a secret
    service_id: str = 'test'       # str or int? Obtained onces for each separate "analysis"
    skyportal_url: str = 'url'
    provenance_origin: str = 'AMPEL' # Will be propagated as part of provenance. Possibly replace with e.g. 'AMPEL@DESY' 

    # Need some mechanism for determining MySourceName
    # name_filter: dict[str, str] = {'ZTF name': 'ZTF', 'TNS ID': 'TNS'}
 
    # Mapping of t2units and contents to parse for submission data. 
    t2_publish: dict[str, Any]

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
        
        Do we more than the the T2DocView here?
        """
        prov = {'origin': self.provenance_origin,
                'stock': t2view.stock, 'link':t2view.link, 'code': t2view.code, 't2_type': t2view.t2_type}
        # More complex - for config we would like the hash rather then full dict
        print('TODO rehash config if needed?')
        prov['config'] = t2view.config
        
        prov['ts'] = t2view.get_time_updated()
        
        return prov
        
    def process(self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None) -> UBson | UnitResult:
        """
        Loop through provided TransientViews, extract data and publish according to the
        configured schema.
        """
        
        for tran_view in gen:

            # Chategorized by stock (internal ampel ID) and channel
            assert tran_view.stock is not None
            stock = tran_view.id
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
                    # If ext_prov is already set, something has gone wrong (should only get here once)
                    if 'external_provenance' in payload['analysis']['results']['data']:
                        print('WARNING+TODO - should exist/log with reasonable message') 

                    # Check that code is good?
                    payload['analysis']['results']['data']['external_provenance'] = ext_prov

                                    
                    # Have some results, check what already is uploaded 
                    print('TODO: Check skyportal for this {object,service}, get our timestap of latest submission')
                
                    # Check whether the latest timestamp of this T2Record is later than what has been published 
                    print('TODO: Published timestamp older than we have, otherwise we would have continued here.')
                    db_ts = t2view.get_time_updated()
                    
                    # Collect other info - plots and results
                    body = t2view.get_payload()
                    if body is None:
                        continue 
                        
                    # Plots
                    print('starting going through plots')
                    if publish_settings['publish_plots'] and (plotlist:=body.get('plots')):
                        payload['analysis']['plots'] = []
                        for plot in plotlist:
                            raw_im = decompress(plot['svg'])
                            payload['analysis']['plots'].append(base64.b64encode(raw_im))

		    # Add analysis information                     
                    for label, path in publish_settings['results_mapping'].items():
                        if result := get_by_path(body, path):
                            payload['analysis']['results'][label] = result
                            
                # Payload for this unit ready.
                if self.debug:
                    print('Payload:')
                    payload['analysis']['plots'][0] = payload['analysis']['plots'][0].decode()[0:20]+'.... cut'
                    print(payload)
                else:
                    # Post payload, will go something like (but through client class?)
                    service_url = f"http://<{self.skyportal_url}>:5000/api/obj/MySourceName/analysis_upload/{self.service_id}"
                    r = requests.post(url, headers={'Authorization':   'token '+self.access_token}, json=payload)
                    print('TODO write something to log')
                     


        return None
