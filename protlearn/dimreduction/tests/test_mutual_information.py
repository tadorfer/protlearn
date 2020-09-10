import pytest
import numpy as np
from ..mutual_information import mutual_information

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_mutual_information():
    "Test MI-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute MI-based feature importances
    X_reduced = mutual_information(X, y, top=10)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        0.        ,  3.93222222, 48.61333333, 86.66666667,  5.        ,
        4.98333333,  4.33666667, -1.30888889,  7.91111111, 19.27888889]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        0.22222222,  4.62222222, 48.45888889, 78.22222222,  5.45888889,
        4.73444444,  4.61333333,  5.52222222,  8.73333333, 20.35311111]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        0.        ,  4.73222222, 43.31555556, 90.33333333,  4.74444444,
        5.09222222,  6.01      , -2.48      ,  8.05555556, 14.03666667]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 10)