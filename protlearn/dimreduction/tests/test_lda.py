import pytest
import numpy as np
from ..lda import lda

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_tree_importance():
    "Test LDA-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute LDA-based feature importances
    X_reduced = lda(X, y)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], -0.57463004)
    np.testing.assert_almost_equal(X_reduced[300,:], 0.12987912)
    np.testing.assert_almost_equal(X_reduced[-1,:], -1.25676961)

    # test array shape
    X_reduced.shape == (700, 1)