import pytest
from ..remove_duplicates import remove_duplicates

def test_remove_duplicates():
    "Test duplicate sequences"

    # define data
    x0 = ['ARKLY', 'LYLPGG', 'EECCKHR']
    x1 = ['ARKLY', 'LYLPGG', 'ARKLY', 'EECCKHR', 'LYLPGG']
    x2 = ['ARKLY', 'LYLPGG', 'ARKLY', 'ARKLY']

    # test for duplicates
    y0 = remove_duplicates(x0)
    y1 = remove_duplicates(x1)
    y2 = remove_duplicates(x2, verbose=2) # checking verbosity

    assert set(y0) == set(x0)
    assert set(y1) == set(x0)
    assert set(y2) == set(['ARKLY', 'LYLPGG'])