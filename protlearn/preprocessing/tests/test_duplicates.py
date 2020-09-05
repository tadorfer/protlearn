import pytest
from ..duplicates import duplicates

def test_duplicates():
    "Test duplicate sequences"

    # define data
    x0 = ['ARKLY', 'LYLPGG', 'EECCKHR']
    x1 = ['ARKLY', 'LYLPGG', 'ARKLY', 'EECCKHR', 'LYLPGG']
    x2 = ['ARKLY', 'LYLPGG', 'ARKLY', 'ARKLY']

    # test for duplicates
    y0 = duplicates(x0)
    y1 = duplicates(x1)
    y2 = duplicates(x2)

    assert set(y0) == set(x0)
    assert set(y1) == set(x0)
    assert set(y2) == set(['ARKLY', 'LYLPGG'])