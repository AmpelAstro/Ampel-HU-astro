
import pytest, logging, re
from unittest.mock import MagicMock
from ampel.contrib.hu.t3.tns import TNSClient, TNSMirrorDB, TNSName

@pytest.mark.skip(reason="requires local TNS mirror")
def test_query_tnsname(t3_transient_views, mocker):
    t3 = TNSMatcher(logging.getLogger())

    assert t3.add(t3_transient_views) is None

    collection = MagicMock()
    mocker.patch('ampel.core.AmpelDB.AmpelDB.get_collection').return_value = collection
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

def test_tnsname():
    name = '2019tdo'
    id_name = TNSName.from_str(name)
    assert id_name.year == 2019
    assert str(id_name) == name
    assert str(TNSName.from_index(int(id_name))) == name

    for i in range(1,4):
        assert TNSName.from_str('2019'+('a'*i)).number == sum(TNSName.BASE**p for p in range(i))
        assert len(str(TNSName(2019,TNSName.BASE**i))) == 4+i
