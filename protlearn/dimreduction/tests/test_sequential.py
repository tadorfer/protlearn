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
    X_fwd = sequential(X, y, estimator=clf)
    X_bwd = sequential(X, y, estimator=clf, direction='backward')

    # test shapes
    X_fwd.shape == (700, 10)
    X_bwd.shape == (700, 10)