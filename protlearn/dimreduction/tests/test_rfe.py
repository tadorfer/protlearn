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
    X = X[:,:50]
    y = np.load(PATH+'features_largeN_labels.npy')

    # perform rfe
    svc = SVC(kernel='linear')
    X_reduced, rank = rfe(X, y, estimator=svc, n_features=10, step=2)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
         0.        , 0.11111111, 0.11111111, 0.        , 0.        ,
         0.        , 0.11111111, 0.71888889, 0.52777778, 0.77]), 
         decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
         0.        , 0.11111111, 0.        , 0.        , 0.11111111,
         0.        , 0.11111111, 0.76333333, 0.56555556, 0.76]), 
         decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
         0.        , 0.        , 0.        , 0.22222222, 0.        ,
         0.        , 0.        , 0.74222222, 0.58666667, 0.66444444]), 
         decimal=3)

    # test shapes
    X_reduced.shape == (700, 10)
    rank.shape == (X.shape[1],)