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
    name="ampel-contrib-hu",
    version="0.7.1-alpha.1",
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
        'ampel-interface>=0.7.1,<0.8',
        'ampel-core[plotting]>=0.7.1,<0.8',
        'ampel-photometry>=0.7.1,<0.8',
        'ampel-alerts>=0.7.1,<0.8',
        'ampel-ztf>=0.7.1,<0.8',
        "catsHTM",
        "extcats",
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
