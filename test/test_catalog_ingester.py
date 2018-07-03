
import pytest, json
import subprocess

from ampel.contrib.hu.TargetCatalogIngester import TargetCatalogIngester
from astropy import units as u
from astropy.time import Time
from pymongo import MongoClient

pytest_plugins = ['ampel.test.fixtures']

@pytest.fixture
def catalog_db(mongod):
	yield mongod
	MongoClient(mongod).drop_database("ToO")

def test_double_init(catalog_db):
	TargetCatalogIngester(catalog_db, "ToO", "testy", 5*u.degree)
	TargetCatalogIngester(catalog_db, "ToO", "testy", 5*u.degree)
	assert MongoClient(catalog_db)["ToO"]["meta"].count({}) == 1

def test_insert(catalog_db):
	jester = TargetCatalogIngester(catalog_db, "ToO", "testy", 5*u.degree)
	t0 = Time.now()
	ra, dec = 100*u.degree, 51*u.degree
	# insert twice to check duplicate behavior
	jester.add_target(ra, dec, 2*u.degree, t0 - 5*u.day, t0 + 25*u.day)
	jester.add_target(ra, dec, 2*u.degree, t0 - 5*u.day, t0 + 25*u.day)

	coll = MongoClient(catalog_db)["ToO"]["testy"]

	count = coll.count({'pos': {'$geoWithin': {"$centerSphere": [[(ra+1*u.deg).to(u.deg).value, (dec).to(u.deg).value], (2*u.deg).to(u.radian).value]}}})
	assert count == 2
