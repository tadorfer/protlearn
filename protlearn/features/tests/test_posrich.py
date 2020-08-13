import pytest
import numpy as np
from features import posrich
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_posrich():
    "Test positional enrichment"

    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'

    # test posrich single position
    posrich_single = posrich(X_list, 2, 'A')
    assert np.array_equal(posrich_single, np.array([1.,0.,1.]))

    # test posrich multiple positions
    posrich_multiple = posrich(X_list, [2, 4], ['A', 'K'])
    assert np.array_equal(posrich_multiple[:,0], np.array([1.,0.,1.]))
    assert np.array_equal(posrich_multiple[:,1], np.array([1.,0.,0.]))
    
    # test ValueError
    with pytest.raises(ValueError):
        posrich_err = posrich(X_err)