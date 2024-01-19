#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-HU-astro/ampel/t0/HealpixTokenConsumer.py
# License           : BSD-3-Clause
# Author            : jn <jno@physik.hu-berlin.de>
# Date              : 28.03.2023
# Last Modified Date: 28.03.2023
# Last Modified By  : jn <jno@physik.hu-berlin.de>

import sys
from signal import SIGINT, SIGTERM, default_int_handler, signal
from typing import Any

import numpy as np
from ampel.alert.AlertConsumer import AlertConsumer
from ampel.alert.AlertConsumerError import AlertConsumerError
from ampel.alert.AlertConsumerMetrics import AlertConsumerMetrics, stat_time
from ampel.contrib.hu.util.AmpelHealpix import AmpelHealpix

# from ampel.util.freeze import recursive_unfreeze
# from ampel.enum.EventCode import EventCode
# from ampel.model.UnitModel import UnitModel
from ampel.core.EventHandler import EventHandler

# from ampel.dev.DevAmpelContext import DevAmpelContext
# from ampel.abstract.AbsAlertSupplier import AbsAlertSupplier
# from ampel.abstract.AbsEventUnit import AbsEventUnit
# from ampel.base.AuxUnitRegister import AuxUnitRegister
# from ampel.alert.FilterBlocksHandler import FilterBlocksHandler
from ampel.ingest.ChainedIngestionHandler import ChainedIngestionHandler
from ampel.log import VERBOSE, AmpelLogger, LogFlag

# from ampel.log.utils import report_exception
from ampel.log.AmpelLoggingError import AmpelLoggingError
from ampel.log.LightLogRecord import LightLogRecord
from ampel.model.ingest.CompilerOptions import CompilerOptions
from ampel.mongo.update.DBUpdatesBuffer import DBUpdatesBuffer

# from ampel.core.AmpelContext import AmpelContext
from ampel.util.mappings import get_by_path, merge_dict
from pymongo.errors import PyMongoError

# from ampel.model.AlertConsumerModel import AlertConsumerModel


class HealpixTokenConsumer(AlertConsumer):
    """
    Complement standard AlertConsumer process by adding Healpix probabilities.

    Normal mode assumes map_name corresponds to resources from HealpixTokenGenerator

    """

    ## Healpix map settings
    map_name: str
    map_hash: None | str = None
    map_url: None | str = None
    map_dir: None | str = None

    # Overload
    def proceed(self, event_hdlr: EventHandler) -> int:
        """
        Process alerts using internal alert_loader/alert_supplier

        :returns: Number of alerts processed
        :raises: LogFlushingError, PyMongoError
        """

        # Setup stats
        #############

        stats = AlertConsumerMetrics(self._fbh.chan_names)

        event_hdlr.set_tier(0)
        run_id = event_hdlr.get_run_id()

        # Setup logging
        ###############

        logger = AmpelLogger.from_profile(
            self.context,
            self.log_profile,
            run_id,
            base_flag=LogFlag.T0 | LogFlag.CORE | self.base_log_flag,
        )

        self.alert_supplier.set_logger(logger)

        if event_hdlr.resources:
            for k, v in event_hdlr.resources.items():
                self.alert_supplier.add_resource(k, v)

        if logger.verbose:
            logger.log(VERBOSE, "Pre-run setup")

        # DBLoggingHandler formats, saves and pushes log records into the DB
        if db_logging_handler := logger.get_db_logging_handler():
            db_logging_handler.auto_flush = False

        # Collects and executes pymongo.operations in collection Ampel_data
        updates_buffer = DBUpdatesBuffer(
            self._ampel_db,
            run_id,
            logger,
            error_callback=self.set_cancel_run,
            catch_signals=False,  # we do it ourself
            max_size=self.updates_buffer_size,
        )

        any_filter = any([fb.filter_model for fb in self._fbh.filter_blocks])

        # Retrieve Healpix map properties
        if map_resource := event_hdlr.resources.get(self.map_name):
            self.map_hash = map_resource.value["hash"]
            self.map_dir = map_resource.value["map_dir"]
        healpix_map = AmpelHealpix(
            map_name=self.map_name, map_url=self.map_url, save_dir=self.map_dir
        )
        map_hash = healpix_map.process_map()
        if not map_hash == self.map_hash:
            raise ValueError("Healpix hash changed - modified map?")

        # Setup ingesters
        ing_hdlr = ChainedIngestionHandler(
            self.context,
            self.shaper,
            self.directives,
            updates_buffer,
            run_id,
            tier=0,
            logger=logger,
            database=self.database,
            trace_id={"alertconsumer": self._trace_id},
            compiler_opts=self.compiler_opts or CompilerOptions(),
        )

        # Loop variables
        iter_max = self.iter_max
        if self.iter_max != self._defaults["iter_max"]:
            logger.debug(f"Using custom iter_max: {self.iter_max}")

        self._cancel_run = 0
        iter_count = 0
        err = 0

        assert self._fbh.chan_names is not None
        reduced_chan_names: str | list[str] = (
            self._fbh.chan_names[0]
            if len(self._fbh.chan_names) == 1
            else self._fbh.chan_names
        )
        fblocks = self._fbh.filter_blocks

        if any_filter:
            filter_results: list[tuple[int, bool | int]] = []
        else:
            filter_results = [(i, True) for i, fb in enumerate(fblocks)]

        # Builds set of stock ids for autocomplete, if needed
        self._fbh.ready(logger, run_id)

        # Shortcuts
        report_filter_error = lambda e, alert, fblock: self._report_ap_error(
            e,
            event_hdlr,
            logger,
            extra={"a": alert.id, "section": "filter", "c": fblock.channel},
        )

        report_ingest_error = lambda e, alert, filter_results: self._report_ap_error(
            e,
            event_hdlr,
            logger,
            extra={
                "a": alert.id,
                "section": "ingest",
                "c": [self.directives[el[0]].channel for el in filter_results],
            },
        )

        # Process alerts
        ################

        # The extra is just a feedback for the console stream handler
        logger.log(self.shout, "Processing alerts", extra={"r": run_id})

        try:
            updates_buffer.start()
            chatty_interrupt = self.chatty_interrupt
            register_signal = self.register_signal

            rejected_count = 0

            # Iterate over alerts
            for alert in self.alert_supplier:
                # Allow execution to complete for this alert (loop exited after ingestion of current alert)
                signal(SIGINT, register_signal)
                signal(SIGTERM, register_signal)

                # Associate upcoming log entries with the current transient id
                stock_id = alert.stock

                if any_filter:
                    filter_results = []

                    # Loop through filter blocks
                    for fblock in fblocks:
                        try:
                            # Apply filter (returns None/False in case of rejection or True/int in case of match)
                            res = fblock.filter(alert)
                            if res[1]:
                                filter_results.append(res)  # type: ignore[arg-type]
                            else:
                                # print("Alert rejected.")
                                rejected_count += 1

                        # Unrecoverable (logging related) errors
                        except (PyMongoError, AmpelLoggingError) as e:
                            print("%s: aborting run() procedure" % e.__class__.__name__)
                            report_filter_error(e, alert, fblock)
                            raise e

                        # Possibly tolerable errors (could be an error from a contributed filter)
                        except Exception as e:
                            if db_logging_handler:
                                fblock.forward(
                                    db_logging_handler,
                                    stock=stock_id,
                                    extra={"a": alert.id},
                                )

                            report_filter_error(e, alert, fblock)

                            if self.raise_exc:
                                raise e
                            else:
                                if self.error_max:
                                    err += 1
                                if err == self.error_max:
                                    logger.error(
                                        "Max number of error reached, breaking alert processing"
                                    )
                                    self.set_cancel_run(
                                        AlertConsumerError.TOO_MANY_ERRORS
                                    )
                else:
                    # if bypassing filters, track passing rates at top level
                    for counter in stats.filter_accepted:
                        counter.inc()

                if filter_results:
                    stats.accepted.inc()

                    # Determine p-value for being associated with Healpix map
                    # Also - should handle missing coord valus gracefully
                    if (
                        len(
                            pos := alert.get_tuples(
                                "ra",
                                "dec",
                                filters=[
                                    {
                                        "attribute": "magpsf",
                                        "operator": "is not",
                                        "value": None,
                                    }
                                ],
                            )
                        )
                        > 0
                    ):
                        ra = np.mean([p[0] for p in pos])
                        dec = np.mean([p[1] for p in pos])
                        alert_cumprob = healpix_map.get_cumprob(ra, dec)
                    else:
                        alert_cumprob = None

                    # Collect information for alert extra
                    alert_extra: dict[str, Any] = {
                        "alert": alert.id,
                        "healpix": {
                            "cumprob": alert_cumprob,
                            "trigger_time": healpix_map.trigger_time,
                            "map_name": self.map_name,
                            "map_hash": self.map_hash,
                            "map_dist": healpix_map.dist,
                            "map_dist_unc": healpix_map.dist_unc,
                            # "map_area": healpix_map.get_maparea()
                        },
                    }

                    try:
                        if self.include_alert_extra_with_keys and alert.extra:
                            for key, path in self.include_alert_extra_with_keys.items():
                                alert_extra[key] = get_by_path(alert.extra, path)
                        with stat_time.labels("ingest").time():
                            ing_hdlr.ingest(
                                alert.datapoints,
                                filter_results,
                                stock_id,
                                alert.tag,
                                alert_extra,
                                alert.extra.get("stock") if alert.extra else None,
                            )
                    except (PyMongoError, AmpelLoggingError) as e:
                        print("%s: abording run() procedure" % e.__class__.__name__)
                        report_ingest_error(e, alert, filter_results)
                        raise e

                    except Exception as e:
                        report_ingest_error(e, alert, filter_results)

                        if self.raise_exc:
                            raise e

                        if self.error_max:
                            err += 1

                        if err == self.error_max:
                            logger.error(
                                "Max number of error reached, breaking alert processing"
                            )
                            self.set_cancel_run(AlertConsumerError.TOO_MANY_ERRORS)

                else:
                    # All channels reject this alert
                    # no log entries goes into the main logs collection sinces those are redirected to Ampel_rej.

                    # So we add a notification manually. For that, we don't use logger
                    # cause rejection messages were alreary logged into the console
                    # by the StreamHandler in channel specific RecordBufferingHandler instances.
                    # So we address directly db_logging_handler, and for that, we create
                    # a LogDocument manually.
                    lr = LightLogRecord(logger.name, LogFlag.INFO | logger.base_flag)
                    lr.stock = stock_id
                    lr.channel = reduced_chan_names  # type: ignore[assignment]
                    lr.extra = {"a": alert.id, "allout": True}
                    if db_logging_handler:
                        db_logging_handler.handle(lr)

                iter_count += 1
                stats.alerts.inc()

                updates_buffer.check_push()
                if db_logging_handler:
                    db_logging_handler.check_flush()

                if iter_count == iter_max:
                    logger.info("Reached max number of iterations")
                    break

                # Exit if so requested (SIGINT, error registered by DBUpdatesBuffer, ...)
                if self._cancel_run > 0:
                    break

                # Restore system default sig handling so that KeyBoardInterrupt
                # can be raised during supplier execution
                signal(SIGINT, chatty_interrupt)
                signal(SIGTERM, chatty_interrupt)

            # print("Rejected count: ", rejected_count)

        # Executed if SIGINT was sent during supplier execution
        except KeyboardInterrupt:
            pass

        except Exception as e:
            event_hdlr.handle_error(e, logger)

            if self.raise_exc:
                raise e

        # Also executed after SIGINT and SIGTERM
        finally:
            updates_buffer.stop()

            if self._cancel_run > 0:
                print("")
                logger.info("Processing interrupted")
            else:
                logger.log(self.shout, "Processing completed")

            try:
                # Flush loggers
                logger.flush()

                # Flush registers and rejected log handlers
                self._fbh.done()

            except Exception as e:
                event_hdlr.handle_error(e, logger)
                if self.raise_exc:
                    raise e

        if self.exit_if_no_alert and iter_count == 0:
            sys.exit(self.exit_if_no_alert)

        # Return number of processed alerts
        return iter_count
