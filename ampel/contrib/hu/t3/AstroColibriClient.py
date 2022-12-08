#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/AstroColibriClient.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                13.11.2022
# Last Modified Date:  13.11.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from typing import Sequence, Dict, Any, Union

import requests
#import json
import backoff
from requests import HTTPError
from requests.auth import HTTPBasicAuth
#from datetime import datetime
#from astropy.time import Time

url = 'https://astro-colibri.herokuapp.com/add_ampel_transient'

class AstroColibriPhot:
    mag: float
    band: str
    datetime: str # e.g. datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")

class AstroColibriPost:
    timestamp: str            # Timestamp of first detection
    trigger_id: str           # AstroColibri identifier -> should be TNS ID
    # ampel_id: str             # Internal AMPEL ID. To be added
    # ivorn: str | None = None  # If existing, this can be added.
    ra: float
    dec: float
    err: float
    type: str = 'ot_sn'       # 'ot_sn'
    observatory: str | None = 'ztf'
    source_name: str          #  SN|AT[space]2022ajn
    ampel_attributes: list[str] = []
    photometry:  dict[str, AstroColibriPhot]
    discovery_name: str      #  ZTFxxx


class AstroColibriClient:
    """
    Initiate a session for communicating with an AstroColibri instance.

    """

    # url (currently test version)
    colibri_url = 'https://astro-colibri.herokuapp.com/add_ampel_transient'


    def __init__(self, colibri_username: str, colibri_password: str, logger):
        self.logger = logger

        self.session = requests.session()
        self.session.auth = HTTPBasicAuth(colibri_username, colibri_password)
        # or ?
        # self.auth=HTTPBasicAuth(colibri_username, colibri_password)


    # robustify post
    @backoff.on_exception(
        backoff.expo,
        requests.ConnectionError,
        max_tries=5,
        factor=10,
        )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: not isinstance(e, HTTPError) or e.response.status_code not in {503, 504, 429, 408},
        max_time=60,
        )
    def firestore_post(self, data: AstroColibriPost )->Dict[Any,Any]:

        print('Post input', data)

        # TODO: not working
        response = self.session.post(self.colibri_url, json=data)
                                #, headers=self.csrfheader)

        #request = requests.post(self.colibri_url, json=data, auth=self.session.auth)

        print('Post response ok ', response.ok)
        print('post response', response.json())


        if response.ok:
            self.logger.debug('AstroColibriClient submit success', extra={"payload": data})
            return {'success':True, **response.json()}

        self.logger.info('AstroColibriClient submit fail', extra={"payload": classification})
        return {'success':False, 'response':response.status_code}
