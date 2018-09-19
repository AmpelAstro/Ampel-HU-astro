
import logging
import pytest
import json
import pkg_resources
from os.path import dirname, join

from ampel.contrib.hu.t2.T2MarshalMonitor import T2MarshalMonitor, ZTFUtils
from ampel.utils.json import object_hook

@pytest.fixture
def lightcurve():
	with open(join(dirname(__file__), 'lightcurve.ZTF18abmjvpb.json')) as f:
		return json.load(f, object_hook=object_hook)

@pytest.fixture
def run_config():
	for resource in pkg_resources.iter_entry_points('ampel.pipeline.t2.configs'):
		if resource.name == 'hu':
			break
	else:
		raise ValueError
	return resource.resolve()()['MARSHALMONITOR_simple']['parameters']

def test_ztf_name(lightcurve):
	ids = lightcurve.get_values('tranId')
	assert len(set(ids)) == 1
	assert ZTFUtils.to_ztf_id(ids[0]) == 'ZTF18abmjvpb'

def test_marshalmonitor(lightcurve, run_config):
	log = logging.getLogger()
	monitor = T2MarshalMonitor(base_config=None, logger=log)
	result = monitor.run(lightcurve, run_config)
	assert result == {'classification': 'SN Ia', 'programs': ['COSM', 'RCF'], 'redshift': '0.066077'}
