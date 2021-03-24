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
        'ampel-interface>=0.7.1-alpha.8,<0.7.2',
        'ampel-core[plotting]>=0.7.1-alpha.7,<0.7.2',
        'ampel-photometry>=0.7.1-alpha.1,<0.7.2',
        'ampel-alerts>=0.7.1-alpha.2,<0.7.2',
        'ampel-ztf>=0.7.1-alpha.10,<0.7.2',
        "catsHTM",
        "extcats",
        # see: https://github.com/sncosmo/sncosmo/issues/291
        "sncosmo==2.2.0",
        "iminuit==1.5.1",
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
