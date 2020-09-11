import pytest
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from ..sequential import sequential

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_sequential():
    "Test sequential feature selection"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    X = X[:,:20]
    y = np.load(PATH+'features_largeN_labels.npy')

    # perform SFS
    clf = RandomForestClassifier(n_estimators=100)
    X_fwd = sequential(X, y, clf)
    X_bwd = sequential(X, y, clf, direction='backward')

    # test array contents
    np.testing.assert_almost_equal(X_fwd[0,:], np.array([
         0.        , 0.11111111, 0.11111111, 0.        , 0.        ,
         0.        , 0.        , 0.11111111, 0.        , 0.]), 
         decimal=3)

    np.testing.assert_almost_equal(X_fwd[-1,:], np.array([
         0.        , 0.        , 0.        , 0.22222222, 0.11111111,
         0.        , 0.        , 0.        , 0.22222222, 0.]), 
         decimal=3)

    np.testing.assert_almost_equal(X_bwd[0,:], np.array([
         0.        , 0.        , 0.11111111, 0.11111111, 0.        ,
         0.        , 0.        , 0.        , 0.        , 0.]), 
         decimal=3)

    np.testing.assert_almost_equal(X_bwd[-1,:], np.array([
         0.        , 0.22222222, 0.        , 0.        , 0.11111111,
         0.22222222, 0.11111111, 0.        , 0.        , 0.]), 
         decimal=3)

    # test shapes
    X_fwd.shape == (700, 10)
    X_bwd.shape == (700, 10)