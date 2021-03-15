import pytest
import numpy as np
from ..ctd import ctd 
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_ctd():
    "Test conjoint triad descriptors"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get ctd
    ctd_arr, desc = ctd(X_list)
    
    # test ctd
    assert ctd_arr.shape == (3, 343)
    assert len(desc) == 343
    assert sum(ctd_arr[0]) == 5
    assert sum(ctd_arr[1]) == 7
    assert sum(ctd_arr[2]) == 6

    # test ValueError
    with pytest.raises(ValueError):
        ctd_error, desc = ctd(X_err)