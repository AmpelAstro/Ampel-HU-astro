
import pytest, logging, re
from unittest.mock import MagicMock
from ampel.contrib.hu.t3.aiotns import TNSMatcher

def test_tnsname(t3_transient_views, mocker):
    t3 = TNSMatcher(logging.getLogger())

    assert t3.add(t3_transient_views) is None

    collection = MagicMock()
    mocker.patch('ampel.db.AmpelDB.AmpelDB.get_collection').return_value = collection
    t3.done()

    assert collection.bulk_write.called_once()
    ops = collection.bulk_write.call_args[0][0]
    assert len(ops) > 0
    pattern = re.compile('^TNS20[0-9]{2}[a-z]{3}$')
    for op in ops:
        assert '$addToSet' in op._doc
        assert 'tranNames' in op._doc['$addToSet']
        assert '$each' in op._doc['$addToSet']['tranNames']
        assert len(op._doc['$addToSet']['tranNames']['$each']) > 0
        for name in op._doc['$addToSet']['tranNames']['$each']:
            assert pattern.match(name)
