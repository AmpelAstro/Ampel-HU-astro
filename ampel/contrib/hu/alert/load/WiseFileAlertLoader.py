#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-alerts/ampel/alert/load/WiseFileAlertLoader.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                30.04.2018
# Last Modified Date:  11.08.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from io import BytesIO
from typing import List
from ampel.abstract.AbsAlertLoader import AbsAlertLoader
import gzip
import json
from pathlib import Path


class WiseFileAlertLoader(AbsAlertLoader[BytesIO]):
    """
    Load alerts from one of more files.
    """

    #: paths to files to load
    file: str

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        if not self.file:
            raise ValueError("Parameter 'files' cannot be empty")

        if self.logger:
            self.logger.info(f"Registering {len(self.file)} file(s) to load")

        if Path(self.file).suffix == ".json":
            self.lc = json.loads(open(self.file, "r").read())
            self.lc_content = iter(self.lc.items())
        elif Path(self.file).suffix == ".gz":
            self.lc = json.loads(gzip.open(self.file, "r").read())
            self.lc_content = iter(self.lc.items())

    def __iter__(self):
        return self

    def __next__(self) -> BytesIO:
        output = json.dumps(next(self.lc_content)).encode("utf-8")
        return BytesIO(output)
