import pytest
from ..remove_unnatural import remove_unnatural

def test_remove_unnatural():
    "Test removal of sequences containing unnatural amino acids"

    # define data
    X = ['ARKLY', 'QERKLI', 'YLLGPGB']

    # test for duplicates
    Y = remove_unnatural(X)

    assert Y == ['ARKLY', 'QERKLI']