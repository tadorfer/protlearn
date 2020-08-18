import pytest
import numpy as np
from features import ctd 
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_paac():
    "Test conjoint triad descriptors"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get ctd
    ctd_list, desc = ctd(X_list)
    
    # test ctd
    assert np.array_equal(ctd_list, np.array([
       [1., 1., 1., 1., 1., 0., 0., 0., 0., 
        0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 1., 1., 1., 1.,
        1., 1., 1., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 
        0., 0., 0., 1., 1., 1., 1., 1., 1.]]))

    # test for longer triads
    assert np.array_equal(ctd('AAAARKLY'), np.array([[2., 1., 1., 1., 1.]]))

    # test ValueError
    with pytest.raises(ValueError):
        ctd_error, desc = ctd(X_err)