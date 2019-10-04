#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/setup.py
# License           : BSD-3-Clause
# Author            : jvs
# Date              : Unspecified
# Last Modified Date: 04.10.2019
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from setuptools import setup, find_packages

setup(
	name='ampel-contrib-hu',
	version='0.7.0',
	# include all packages under 'ampel'
	packages=find_packages('ampel'),
	# if any package contains *.json files, include them
	package_data = {'': ['*.json']},
	entry_points = {
		'ampel_target_sources' : [
			'TargetSourceListener = ampel.contrib.hu.TargetSourceListener:TargetSourceListener',
			'TargetSourceListenerSlack = ampel.contrib.hu.TargetSourceListenerSlack:TargetSourceListenerSlack',
		],
		'ampel_resources' : [
			'extcats = ampel.contrib.hu.resources:extcatsURI',
			'catsHTM = ampel.contrib.hu.resources:catsHTMURI',
			'desycloud = ampel.contrib.hu.resources:desyCloudURI',
		],
		'console_scripts' : [
			'catshtmd = ampel.contrib.hu.catshtm_server:run'
		]
	}
)
