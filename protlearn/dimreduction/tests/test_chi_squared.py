import pytest
import numpy as np
from ..chi_squared import chi_squared

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_chi_squared():
    "Test chi-squared-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute chi-squared-based feature importances
    X_reduced = chi_squared(X, y, top=10)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        0.        , 0.        , 0.        , 0.        , 0.19548082,
        0.49728261, 0.1206297 , 0.23030149, 0.10641248, 0.44082605]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        0.        , 0.        , 0.16666667, 0.        , 0.21650026,
        0.61277174, 0.17247985, 0.02685584, 0.15840555, 0.52223987]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        0.5       , 0.2       , 0.        , 0.        , 0.20651603,
        0.57540761, 0.43422943, 0.48441855, 0.36048527, 0.1457506]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 10)