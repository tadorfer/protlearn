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

    # test shapes
    X_reduced.shape == (700, 10)
    rank.shape == (X.shape[1],)