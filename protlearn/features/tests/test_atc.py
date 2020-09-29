import pytest
import numpy as np
from ..atc import atc
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_atc():
    "Test sequence compositions"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test relative composition
    atc_rel, bonds = atc(X_list)

    # test array contents
    np.testing.assert_almost_equal(atc_rel[0], np.array([
        0.27083333, 0.54861111, 0.07638889, 0.10416667, 0.]))
    np.testing.assert_almost_equal(atc_rel[1], np.array([
        0.2585034 , 0.52380952, 0.06122449, 0.14965986, 0.00680272]))
    np.testing.assert_almost_equal(atc_rel[2], np.array([
        0.234375 , 0.5      , 0.09375  , 0.1640625, 0.0078125]))

    np.testing.assert_almost_equal(bonds[0], np.array([138., 127.,  11.]))
    np.testing.assert_almost_equal(bonds[1], np.array([140., 129.,  11.]))
    np.testing.assert_almost_equal(bonds[2], np.array([120., 108.,  12.])) 
    
    # test if frequencies add to 1
    for i in range(len(X_list)):
        assert round(atc_rel[0].sum()) == 1
    
    # test absolute composition
    atc_abs, bonds = atc(X_list, method='absolute')

    # test array contents
    assert np.array_equal(atc_abs[0], np.array([
        39., 79., 11., 15.,  0.]))
    assert np.array_equal(atc_abs[1], np.array([
        38., 77.,  9., 22.,  1.]))
    assert np.array_equal(atc_abs[2], np.array([
        30., 64., 12., 21.,  1.]))

    # test ValueError
    with pytest.raises(ValueError):
        atc_err, bond_err = atc(X_err)