
import pytest
from ampel.contrib.hu.t3.aiotns import TNSAnnotate
pytest_plugins = ['ampel.test.fixtures']

def test_tnsname(t3_transient_views, mocker):
    t3 = TNSAnnotate(None)

    journal = t3.add(t3_transient_views)
    assert len(journal) > 0
    for entry in journal:
        assert 'tnsNames' in entry.content
        assert len(entry.content['tnsNames']) > 0
        for name in entry.content['tnsNames']:
            assert name.startswith('AT') or name.startswith('SN')
    
    t3.done()

