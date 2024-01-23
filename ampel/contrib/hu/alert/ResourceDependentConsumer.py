#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-HU-astro/ampel/alert/ResourceDependentConsumer.py
# License           : BSD-3-Clause
# Author            : jno
# Date              : 28.03.2023
# Last Modified Date: 28.03.2023
# Last Modified By  : jno

import json
import re

from ampel.alert.AlertConsumer import AlertConsumer
from ampel.core.EventHandler import EventHandler
from ampel.ingest.ChainedIngestionHandler import ChainedIngestionHandler
from ampel.log import AmpelLogger
from ampel.model.ingest.CompilerOptions import CompilerOptions
from ampel.model.ingest.IngestDirective import IngestDirective
from ampel.mongo.update.DBUpdatesBuffer import DBUpdatesBuffer


class DynamicShaperConsumer(AlertConsumer):
    """

    Extension of standard AlertConsumer where the configuration of the shaper can be updated based
    (dynamic) resources available to the EventHandler.

    Use case is when a config parameter is created dynamically in the process of running the job
    and should be added to the alert information stored into the db
    Since updated values are present during configuration building, these will still be part of the
    provenance/hash sequence.

    :param directives_map:

    For each key of directives_map with a value which exists in the resources,
    which corresponds to a resource entry, any occurance of the value
    in the directies will be replaced with the resource value.


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

    directives_map: dict[str, str]

    # Overload
    def get_ingestion_handler(
        self,
        event_hdlr: EventHandler,
        updates_buffer: DBUpdatesBuffer,
        logger: AmpelLogger,
    ) -> ChainedIngestionHandler:
        # Replacement pattern for provided resources
        pattern = re.compile(
            "|".join(
                [
                    k
                    for k, v in self.directives_map.items()
                    if v in self.alert_supplier.resources
                ]
            )
        )
        print("pattern", pattern)
        print("directives pre", self.directives)
        directives = [
            IngestDirective(
                **json.loads(
                    pattern.sub(
                        lambda m: self.alert_supplier.resources[
                            self.directives_map[m.group(0)]
                        ].value,
                        json.dumps(el.dict()),
                    )
                )
            )
            for el in self.directives
        ]
        print("directives post", directives)

        return ChainedIngestionHandler(
            self.context,
            self.shaper,
            directives,
            updates_buffer,
            event_hdlr.get_run_id(),
            tier=0,
            logger=logger,
            database=self.database,
            trace_id={"alertconsumer": self._trace_id},
            compiler_opts=self.compiler_opts or CompilerOptions(),
        )
