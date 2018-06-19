
import pytest, json
import subprocess

from ampel.contrib.hu.TargetCatalogIngester import TargetCatalogIngester
from astropy import units as u
from astropy.time import Time

def docker_service(image, protocol, port):
	container = None
	try:
		container = subprocess.check_output(['docker', 'run', '--rm', '-d', '-P', image]).strip()
		ports = json.loads(subprocess.check_output(['docker', 'container', 'inspect', '-f', '{{json .NetworkSettings.Ports}}', container]))
		yield '{}://localhost:{}'.format(protocol, ports['{}/tcp'.format(port)][0]['HostPort'])
	except FileNotFoundError:
		return pytest.skip("Docker fixture requires Docker")
	finally:
		if container is not None:
			subprocess.check_call(['docker', 'container', 'stop', container],
			    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

@pytest.fixture(scope="session")
def mongod():
	yield from docker_service('mongo:3.6', 'mongodb', 27017)

def test_insert(mongod):
	jester = TargetCatalogIngester(mongod, "ToO", "testy", 5*u.degree)
	t0 = Time.now()
	ra, dec = 100*u.degree, 51*u.degree
	# insert twice to check duplicate behavior
	jester.add_target(ra, dec, 2*u.degree, t0 - 5*u.day, t0 + 25*u.day)
	jester.add_target(ra, dec, 2*u.degree, t0 - 5*u.day, t0 + 25*u.day)
	
	from pymongo import MongoClient
	mc = MongoClient(mongod)
	coll = mc["ToO"]["testy"]
	
	count = coll.count({'pos': {'$geoWithin': {"$centerSphere": [[(ra+1*u.deg).to(u.deg).value, (dec).to(u.deg).value], (2*u.deg).to(u.radian).value]}}})
	assert count == 2
