import pytest
import numpy as np
from ..geary import geary
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_geary():
    "Test Geary's C"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test Geary's C
    g = geary(X_list)

    # test array contents
    np.testing.assert_almost_equal(g[0], np.array([
        0.39007928, 0.7152915 , 0.65898912, 0.84985622, 0.64272388,
        0.4059528 , 0.23863636, 0.17678866]))
    np.testing.assert_almost_equal(g[1], np.array([
        1.19694774, 0.61612546, 0.75825192, 0.60961432, 0.65681427,
        0.65348704, 0.39620939, 1.39464463]))
    np.testing.assert_almost_equal(g[2], np.array([
        1.1618387 , 0.86868123, 0.74366823, 0.5402789 , 0.79049819,
        0.76883367, 0.69269949, 1.37886563]))

    # test ValueError (alphabetical)
    with pytest.raises(ValueError):
        m_err = geary(X_err)

    # test ValueError (maximum lag)
    with pytest.raises(ValueError):
        m_err = geary(X_list, d= 31) 

    # test ValueError (d >= min_len)
    with pytest.raises(ValueError):
        m_err = geary(X_err, d=7)