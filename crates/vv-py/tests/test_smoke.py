import vizviz


def test_version_is_exposed():
    assert vizviz.__version__


def test_element_symbol():
    assert vizviz.element_symbol(6) == "C"
    assert vizviz.element_symbol(26) == "Fe"
    assert vizviz.element_symbol(0) is None
