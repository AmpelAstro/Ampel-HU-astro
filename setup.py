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
    version="0.7.0",
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
        "ampel-interface",
        "ampel-core",
        "ampel-photometry",
        "ampel-alerts",
        "ampel-ztf",
        "catsHTM",
        "extcats",
        "zerorpc",
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
    ],
    entry_points={
        "ampel_target_sources": [
            "TargetSourceListener = ampel.contrib.hu.TargetSourceListener:TargetSourceListener",
            "TargetSourceListenerSlack = ampel.contrib.hu.TargetSourceListenerSlack:TargetSourceListenerSlack",
        ],
        "ampel_resources": [
            "extcats = ampel.contrib.hu.resources:extcatsURI",
            "catsHTM = ampel.contrib.hu.resources:catsHTMURI",
            "desycloud = ampel.contrib.hu.resources:desyCloudURI",
        ],
        "console_scripts": ["catshtmd = ampel.contrib.hu.catshtm_server:run"],
    },
)
