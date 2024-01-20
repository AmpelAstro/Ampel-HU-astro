from ampel.contrib.hu.t3.tns.TNSName import TNSName


def test_tnsname():
    name = "2019tdo"
    id_name = TNSName.from_str(name)
    assert id_name.year == 2019
    assert str(id_name) == name
    assert str(TNSName.from_index(int(id_name))) == name

    for i in range(1, 4):
        assert TNSName.from_str("2019" + ("a" * i)).number == sum(
            TNSName.BASE**p for p in range(i)
        )
        assert len(str(TNSName(2019, TNSName.BASE**i))) == 4 + i
