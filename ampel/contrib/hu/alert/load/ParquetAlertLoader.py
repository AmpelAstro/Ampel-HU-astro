#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/alert/load/ParquetAlertLoader.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen
# Date:                23.02.2026
# Last Modified Date:  02.04.2026
# Last Modified By:    Felix Fischer <felix.martin.fischer@desy.de>

import re
from collections.abc import Iterator
from typing import Any

from ampel.abstract.AbsAlertLoader import AbsAlertLoader


def _get_pyarrow_dataset() -> Any:
    try:
        import pyarrow.dataset as ds  # type: ignore[import-not-found,import-untyped]
    except ImportError as exc:
        raise ImportError(
            "ParquetAlertLoader requires the optional dependency 'pyarrow'. "
            "Please install it to use this loader."
        ) from exc
    return ds


def _parse_condition(condition: str | None) -> Any | None:
    """
    Parse a condition string into an expression usable by pyarrow dataset.
    This is used to parse the filtering conditions. The filter uses the same
    semantics as LSSTArchiveAlertLoader.
    """
    if not condition:
        return None

    ds = _get_pyarrow_dataset()
    parts = re.split(r"\n|\band\b", condition, flags=re.IGNORECASE)

    filt = None

    for part in parts:
        clause = part.strip()
        if not clause:
            continue

        m = re.match(
            r"([A-Za-z0-9_.]+)\s*(>=|<=|!=|==|=|>|<)\s*(.+)",
            clause,
        )
        if not m:
            raise ValueError(f"Invalid condition: {clause}")

        key, op, raw_value = m.groups()
        raw_value = raw_value.strip()

        value: None | bool | int | float | str

        if raw_value.lower() == "null":
            value = None
        elif raw_value.lower() == "true":
            value = True
        elif raw_value.lower() == "false":
            value = False
        else:
            try:
                value = float(raw_value) if "." in raw_value else int(raw_value)
            except ValueError:
                value = raw_value.strip('"').strip("'")

        field = ds.field(tuple(key.split(".")))

        if op in ("=", "=="):
            expr = field == value
        elif op == "!=":
            expr = field != value
        elif op == ">":
            expr = field > value
        elif op == ">=":
            expr = field >= value
        elif op == "<":
            expr = field < value
        elif op == "<=":
            expr = field <= value
        else:
            raise ValueError(f"Unsupported operator: {op}")

        filt = expr if filt is None else (filt & expr)

    return filt


class ParquetAlertLoader(AbsAlertLoader[dict]):
    """
    Alert loader for parquet files. This is used for testing and development purposes.
    """

    path: str
    condition: str | None = None

    def __init__(self, **kwargs) -> None:
        """
        Initialize the ParquetAlertLoader.

        Note: We now use pyarrow.dataset.scanner(filter=...) instead of ParquetFile,
        to enable Arrow-level filtering.
        """
        ds = _get_pyarrow_dataset()

        super().__init__(**kwargs)

        dataset = ds.dataset(self.path, format="parquet")
        filt = _parse_condition(self.condition)

        scanner = dataset.scanner(
            filter=filt,
            batch_size=1,
        )

        self._batches = iter(scanner.to_batches())
        self._batch: Iterator[dict] | None = None

    def _next_batch(self) -> Iterator[dict]:
        """
        Return an iterator over the next batch of rows from the parquet file.

        Because we use batch_size=1, Arrow may yield empty batches after filtering, which
        falsely signals end-of-data. We must skip empty batches explicitly to avoid premature
        StopIteration.
        """
        while True:
            batch = next(self._batches)
            rows = batch.to_pylist()
            if rows:
                return iter(rows)

    def __next__(self) -> dict:
        """
        Return the next alert as a dictionary.

        This method is used to implement the iterator protocol. It is called automatically
        when the iterator is advanced to the next item. If the end of the data is reached, it
        raises a StopIteration exception.
        """
        if self._batch is None:
            self._batch = self._next_batch()
        try:
            return next(self._batch)
        except StopIteration:
            self._batch = self._next_batch()
            return next(self._batch)