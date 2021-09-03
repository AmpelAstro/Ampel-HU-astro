#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : setup.py
# License           : BSD-3-Clause
# Author            : jvs
# Date              : Unspecified
# Last Modified Date: 06.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from setuptools import find_namespace_packages, setup

setup(
    name="ampel-hu-astro",
    version="0.8.0-alpha.0",
    packages=find_namespace_packages(),
    package_data={
        "": ["*.json", "py.typed"],  # include any package containing *.json files
        "conf": [
            "*.json",
            "**/*.json",
            "**/**/*.json",
            "*.yaml",
            "**/*.yaml",
            "**/**/*.yaml",
            "*.yml",
            "**/*.yml",
            "**/**/*.yml",
        ],
    },
    install_requires=[
        'ampel-interface>=0.8.0a0,<0.9',
        'ampel-core>=0.8.0a6,<0.9',
        'ampel-plots>=0.8.0a0,<0.9',
        'ampel-photometry>=0.8.0a0,<0.9',
        'ampel-alerts>=0.8.0a0,<0.9',
        'ampel-ztf>=0.8.0a0,<0.9',
        "sncosmo",
        "iminuit",
        "sfdmap",
        "astropy",
        "numpy",
        "scipy",
        "beautifulsoup4",
        "backoff",
        "requests",
        "pymage @ https://github.com/MickaelRigault/pymage/archive/v1.0.tar.gz",
        # pymage secretly depends on pandas
        "pandas",
    ],
)
