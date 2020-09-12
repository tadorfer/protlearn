import pytest
import numpy as np
from ..univariate_filter import univariate_filter

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_univariate_filter():
    "Test filter-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute filter-based feature importances
    X_f = univariate_filter(X, y, top=10)
    X_chi2 = univariate_filter(X, y, method='chi2', top=10)
    X_mi = univariate_filter(X, y, method='mutual_info', top=10)

    # test f_test
    np.testing.assert_almost_equal(X_f[0,:], np.array([
        0.        ,  0.71888889, -0.32333333,  4.33666667, -1.198     ,
       13.93111111, -0.13444444,  1.07      , -0.07111111,  4.43333333]), decimal=3)

    np.testing.assert_almost_equal(X_f[300,:], np.array([
        0.11111111,  0.76333333, -0.44888889,  4.61333333, -1.13922222,
       13.48111111,  0.02666667,  1.37888889, -0.16033333,  4.58888889]), decimal=3)

    np.testing.assert_almost_equal(X_f[-1,:], np.array([
        0.00000000e+00,  7.42222222e-01,  4.44444444e-03,  6.01000000e+00,
       -1.21611111e+00,  1.43533333e+01,  1.11111111e-03,  2.09777778e+00,
        4.03333333e-02,  6.12222222e+00]), decimal=3)
    
    # test chi2
    np.testing.assert_almost_equal(X_chi2[0,:], np.array([
        0.        , 0.        , 0.        , 0.        , 0.19548082,
        0.49728261, 0.1206297 , 0.23030149, 0.10641248, 0.44082605]), decimal=3)

    np.testing.assert_almost_equal(X_chi2[300,:], np.array([
        0.        , 0.        , 0.16666667, 0.        , 0.21650026,
        0.61277174, 0.17247985, 0.02685584, 0.15840555, 0.52223987]), decimal=3)

    np.testing.assert_almost_equal(X_chi2[-1,:], np.array([
        0.5       , 0.2       , 0.        , 0.        , 0.20651603,
        0.57540761, 0.43422943, 0.48441855, 0.36048527, 0.1457506]), decimal=3)

    # test array shape
    X_f.shape == (700, 10)
    X_chi2.shape == (700, 10)
    X_mi.shape == (700, 10)