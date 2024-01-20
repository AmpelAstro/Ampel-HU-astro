from ampel.contrib.hu.utils import flatten


def test_flatten():
    assert flatten(1, 2, ["b", "a", ["c", "d"]], 3) == [1, 2, "b", "a", "c", "d", 3]
