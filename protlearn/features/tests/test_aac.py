import pytest
import numpy as np
from ..aac import aac
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_aac():
    "Test sequence compositions"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test relative composition
    aac_rel, aa = aac(X_list, remove_zero_cols=True)

    # test array contents
    np.testing.assert_almost_equal(aac_rel[0], np.array([
        0.28571429, 0.        , 0.        , 0.        , 0.        ,
        0.14285714, 0.28571429, 0.        , 0.        , 0.14285714,
        0.14285714]))
    np.testing.assert_almost_equal(aac_rel[1], np.array([
        0.        , 0.11111111, 0.11111111, 0.11111111, 0.22222222,
        0.        , 0.22222222, 0.        , 0.22222222, 0.        ,
        0.]))
    np.testing.assert_almost_equal(aac_rel[2], np.array([
        0.375     , 0.125     , 0.25      , 0.        , 0.        ,
        0.        , 0.        , 0.125     , 0.        , 0.125     ,
        0.]))
    
    # test if frequencies add to 1
    for i in range(len(X_list)):
        assert round(aac_rel[0,:].sum()) == 1
    
    # test absolute composition
    aac_abs, aa = aac(X_list, method='absolute', remove_zero_cols=True)

    # test array contents
    assert np.array_equal(aac_abs[0], np.array([
        2., 0., 0., 0., 0., 1., 2., 0., 0., 1., 1.]))
    assert np.array_equal(aac_abs[1], np.array([
        0., 1., 1., 1., 2., 0., 2., 0., 2., 0., 0.]))
    assert np.array_equal(aac_abs[2], np.array([
        3., 1., 2., 0., 0., 0., 0., 1., 0., 1., 0.]))
    
    # test if frequences == sequence length
    all_lengths = [7,9,8]
    for i in range(len(X_list)):
        assert aac_abs[i,:].sum() == all_lengths[i]

    # test ValueError
    with pytest.raises(ValueError):
        acc_err, aa = aac(X_err)