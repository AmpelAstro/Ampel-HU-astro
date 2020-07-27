import pytest
from os.path import dirname, join
from ampel.util.json import load
import logging


import ampel.contrib.hu.t3.TNSTalker as TNSTalker


# Load mixed set of transients for check
@pytest.fixture
def t3_TNStvs_mixed():
    with open(join(dirname(__file__),'t3TransientTalker_testTV.json')) as f:    
#    with open('t3TransientTalker_testTV.json') as f:
        views = [v for v in load(f)]
    return views


@pytest.fixture
def t3_TNStvs_good():
    with open(join(dirname(__file__),'t3TransientTalker_testsubmitTV.json')) as f:    
#    with open('t3TransientTalker_testsubmitTV.json') as f:
        views = [v for v in load(f)]
    return views



@pytest.fixture
def testrunconfig():
    return TNSTalker.TNSTalker.RunConfig(
        tns_api_key = "a3f9bcbbe6a26a1ae97a0eeefe465932e68cba83",
        submit_tns = False,
        sandbox = False,
        needed_catalogs = []
    )


def test_run_t3_TnsTalker_selection(t3_TNStvs_mixed, testrunconfig):
    '''
    Check whether correct TVs are accepted based on transient information
    '''
    unit_test = TNSTalker.TNSTalker(logger=logging.getLogger(),run_config=testrunconfig)
    out = [unit_test.accept_tview(tv) for tv in t3_TNStvs_mixed]
    assert len(out) == 30
    assert sum(out) == 5



def test_run_t3_TnsTalker_tnsnamefound(t3_TNStvs_good, testrunconfig):
    '''
    Check whether it find all of these to be submitted. Note - can connect to remove TNS DB!
    '''
    unit_test = TNSTalker.TNSTalker(logger=logging.getLogger(),run_config=testrunconfig)
    for tv in t3_TNStvs_good:
        tns_name, tns_internals, jup = unit_test.find_tns_name(tv)
        assert tns_name is not None



