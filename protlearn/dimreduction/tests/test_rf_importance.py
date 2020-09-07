import pytest
import numpy as np
from ..rf_importance import rf_importance

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_rf_importance():
    "Test RF-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')

    # compute RF-based feature importances
    X_reduced = rf_importance(X, y, top=5)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        0.37577778, 0.28177778, 8.36255556, 2.22222222, 6.64444444]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        0.77022222,   0.31077778,   8.37488889, -43., 5.97777778]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        0.60611111,  0.365,  8.36233333, -0.44444444,  6.73333333]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 5)