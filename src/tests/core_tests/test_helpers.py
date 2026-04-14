from amuse.support.helpers import rename_fn_par

import pytest


@pytest.mark.filterwarnings("ignore: old is deprecated.*")
def test_rename_fn_par():
    assert rename_fn_par("new", None, "old", None, 1) == 1
    assert rename_fn_par("new", None, "old", 1, 2) == 1
    assert rename_fn_par("new", 1, "old", None, 2) == 1
    assert rename_fn_par("new", 1, "old", 1, 2) == 1

    with pytest.raises(ValueError):
        rename_fn_par("new", 1, "old", 2, 3)


class _RenamedMethod:
    def __init__(self, new = None, a = 1, b = "test", old = None):
        self.y = rename_fn_par("new", new, "old", old, 42)


@pytest.mark.filterwarnings("ignore: old is deprecated.*")
def test_rename_fn_par_usage():
    rm = _RenamedMethod()
    assert rm.y == 42

    rm = _RenamedMethod(1)
    assert rm.y == 1

    rm = _RenamedMethod(new=1)
    assert rm.y == 1

    rm = _RenamedMethod(old=1)
    assert rm.y == 1

    rm = _RenamedMethod(1, 1, "test", 1)
    assert rm.y == 1

    with pytest.raises(ValueError):
        _RenamedMethod(1, 2, "test", 3)
