import pytest
import numpy as np
from ..rf_importance import rf_importance

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_rf_importance():
    "Test RF-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute RF-based feature importances
    X_reduced = rf_importance(X, y, top=5)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        1.38777778,  4.33666667,  1.12444444, 43.64444444,  0.90533333]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        1.11666667,  4.61333333,  1.17111111, 31.83333333,  0.90266667]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        1.25888889,  6.01      ,  1.08      , 22.86666667,  0.92633333]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 5)