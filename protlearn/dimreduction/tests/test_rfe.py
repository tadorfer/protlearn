import pytest
import numpy as np
from sklearn.svm import SVC
from ..rfe import rfe

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_rfe():
    "Test RFE-based feature selection"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # perform rfe
    svc = SVC(kernel='linear')
    X_reduced, rank = rfe(X, y, svc, n_features=10, step=5)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
         2.25555556, 12.2       , 86.66666667,  4.90111111,  4.98333333,
         19.40444444, -0.27777778,  1.38888889, 10.12222222, 24.66666667]), 
         decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
         3.06666667, 11.21111111, 78.22222222,  4.65444444,  4.73444444,
         18.76777778,  0.67777778,  0.32222222,  8.12222222, 24.44444444]), 
         decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
         1.36666667e+00,  1.21111111e+01,  9.03333333e+01,  5.45666667e+00,
         5.09222222e+00,  2.15166667e+01, -5.77777778e-01, -8.88888889e-02,
         9.76666667e+00,  1.74444444e+01]), 
         decimal=3)

    # test shapes
    X_reduced.shape == (700, 10)
    rank.shape == (X.shape[1],)