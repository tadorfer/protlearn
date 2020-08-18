import pytest
import numpy as np
from features import cksaap
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_paac():
    "Test k-spaced amino acid pair composition"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get paac
    cksaap_list, desc = cksaap(X_list)
    
    # test paac
    assert np.array_equal(cksaap_list, np.array([
       [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0],
       [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
       [1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1]]))

    # test ValueError
    with pytest.raises(ValueError):
        cksaap_error, desc = cksaap(X_err)