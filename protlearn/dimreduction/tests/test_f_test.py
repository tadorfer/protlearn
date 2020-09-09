import pytest
import numpy as np
from ..f_test import f_test

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_f_test():
    "Test f-test-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute f-test-based feature importances
    X_reduced = f_test(X, y, top=10)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
        0.        ,  0.71888889, -0.32333333,  4.33666667, -1.198     ,
       13.93111111, -0.13444444,  1.07      , -0.07111111,  4.43333333]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
        0.11111111,  0.76333333, -0.44888889,  4.61333333, -1.13922222,
       13.48111111,  0.02666667,  1.37888889, -0.16033333,  4.58888889]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
        0.00000000e+00,  7.42222222e-01,  4.44444444e-03,  6.01000000e+00,
       -1.21611111e+00,  1.43533333e+01,  1.11111111e-03,  2.09777778e+00,
        4.03333333e-02,  6.12222222e+00]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 10)