import pytest
import numpy as np
from ..xgb_importance import xgb_importance

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_xgb_importance():
    "Test XGB-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute XGB-based feature importances
    X_reduced = xgb_importance(X, y, top=5)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        -0.01,  1.14444444,  4.73777778,  0.16666667, -0.67777778]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        0.82777778,  1.03333333,  6.01888889,  0.07333333, -0.57888889]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        -0.52222222,  0.95555556,  6.02444444,  0.82555556, -0.97555556]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 5)