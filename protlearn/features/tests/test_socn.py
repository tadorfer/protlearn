import pytest
import numpy as np
from ..socn import socn
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_socn():
    "Test sequence-order-coupling number"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get socn
    socn_sw, socn_g = socn(X_list, d=5)
    
    # test socn
    np.testing.assert_almost_equal(socn_sw, 
    np.array([[1.786946, 3.375387, 2.740649, 1.31594 , 0.329672],
              [3.148735, 2.022538, 0.876211, 2.230154, 0.900203],
              [1.843719, 2.377882, 2.331234, 1.830166, 2.160547]]),
        decimal=3)

    np.testing.assert_almost_equal(socn_g, 
    np.array([[ 21741.,  42454.,  45633.,  32164.,  18432.],
              [117964.,  95881.,  89542.,  94674.,  54093.],
              [67902.,  87135.,  79295.,  77173.,  64152.]]),
        decimal=3)

    # test ValueError (alphabetical)
    with pytest.raises(ValueError):
        socn_error_sw, socn_error_g = socn(X_err, d=5)

    # test ValueError (exceeding maximum lag)
    with pytest.raises(ValueError):
        socn_error_sw, socn_error_g = socn(X_err, d=15)