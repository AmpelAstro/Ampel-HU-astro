#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-alerts/ampel/alert/load/WiseFileAlertLoader.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                30.04.2018
# Last Modified Date:  11.08.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

import gzip
import json
from io import BytesIO
from pathlib import Path

from ampel.abstract.AbsAlertLoader import AbsAlertLoader


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
            with open(self.file) as f:
                self.lc = json.load(f)
            self.lc_content = iter(self.lc.items())
        elif Path(self.file).suffix == ".gz":
            with gzip.open(self.file) as f:
                self.lc = json.load(f)
            self.lc_content = iter(self.lc.items())

    def __iter__(self):
        return self

    def __next__(self) -> BytesIO:
        output = json.dumps(next(self.lc_content)).encode("utf-8")
        return BytesIO(output)
