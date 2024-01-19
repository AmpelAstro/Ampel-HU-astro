#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-HU-astro/ampel/alert/DynamicShaperAlertConsumer.py
# License           : BSD-3-Clause
# Author            : jno
# Date              : 28.03.2023
# Last Modified Date: 28.03.2023
# Last Modified By  : jno

import re
import json

from ampel.log import AmpelLogger
from ampel.alert.AlertConsumer import AlertConsumer
from ampel.model.UnitModel import UnitModel
from ampel.model.ingest.IngestDirective import IngestDirective
from ampel.ingest.ChainedIngestionHandler import ChainedIngestionHandler
from ampel.mongo.update.DBUpdatesBuffer import DBUpdatesBuffer
from ampel.model.ingest.CompilerOptions import CompilerOptions


class DynamicShaperAlertConsumer(AlertConsumer):
    """

    Extension of standard AlertConsumer where the configuration of the shaper can be updated based
    (dynamic) resources available to the EventHandler.

    Use case is when a config parameter is created dynamically in the process of running the job
    and should be added to the alert information stored into the db

    :param shaper_map:

    Transfer values in the shaper_map to a corresponding entry in the alert shaper config.

    Ok, new take on this. Instead of changing the directives (which are already hashed and set), we wish to
    add a new datapoint to represent the map. To do so we minimally have to change the ZiDataPointShaper since
    this only accepts ztf like datapoints (detections or upper limits). We here wish to have a datapoint which is
    shared by all stocks (if possible). So maybe this is where all changes should go?

    The shaper is initaialized by the ChainedIngestionHandler, so should be fine to change config here.
    Steps:
    1. Create version of ZiDataPointShaper which adds GW datapoint.
    2. Change name of this class to something shaper related, then change logic below to instead create
    appropriate config based on the resources.

    """

    shaper_map: dict[str, str]

    # Overload
    def get_ingestion_handler(
        self, run_id: int, updates_buffer: DBUpdatesBuffer, logger: AmpelLogger
    ) -> ChainedIngestionHandler:
        # Update shaper
        # print('config first', self.shaper)
        newconfig = {
            config_key: self.alert_supplier.resources[resource_name].value
            for config_key, resource_name in self.shaper_map.items()
            if resource_name in self.alert_supplier.resources
        }
        # print('config pure', newconfig)
        if isinstance(self.shaper.config, dict):
            newconfig = self.shaper.config | newconfig
        shaper = UnitModel(
            unit=self.shaper.unit,
            config=newconfig,
            secrets=self.shaper.secrets,
            override=self.shaper.override,
        )
        # print('config after', shaper)

        return ChainedIngestionHandler(
            self.context,
            shaper,
            self.directives,
            updates_buffer,
            run_id,
            tier=0,
            logger=logger,
            database=self.database,
            trace_id={"alertconsumer": self._trace_id},
            compiler_opts=self.compiler_opts or CompilerOptions(),
        )
