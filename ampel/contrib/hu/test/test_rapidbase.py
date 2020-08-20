import pytest
from os.path import dirname, join
from ampel.util.json import load
import logging


import ampel.contrib.hu.t3.RapidBase as RapidBase


# Load mixed set of transients for check
@pytest.fixture
def t3_rapid_tvs():
    with open(join(dirname(__file__),'t3RapidBase_testsubmitTV.json')) as f:    
        views = [v for v in load(f)]
    return views


@pytest.fixture
def testrunconfig():
    return RapidBase.RapidBase.RunConfig(
        do_react = False,
        do_testreact = False, 
        redshift_catalogs = ('NEDz', 'NEDz_extcats', 'SDSS_spec'), 
        min_peak_mag = 20.0, 
        max_age = 20,
        max_redshift = 0.8, 
        min_ndet = 1, 
        min_ndet_postul = 1, 
        max_absmag = -10,
        min_absmag = -20,
        maglim_maxago = 8
    )




def test_run_t3_RapidBase_selected(t3_rapid_tvs, testrunconfig):
    '''
    Check rapid selection info output
    '''
    unit_test = RapidBase.RapidBase(logger=logging.getLogger(),run_config=testrunconfig)
    selected = [tv for tv in t3_rapid_tvs if unit_test.accept_tview(tv)]
    assert len(selected) == 3



