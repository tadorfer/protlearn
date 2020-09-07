import pytest
import numpy as np
from ..xgb_importance import xgb_importance

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_xgb_importance():
    "Test XGB-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')

    # compute XGB-based feature importances
    X_reduced = rf_importance(X, y, top=5)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        0.11111111, 0.        , 0.        , 0.        , 0.]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        0.11111111, 0.        , 0.        , 0.        , 0.]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        0.        , 0.        , 0.22222222, 0.        , 0.11111111]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 5)