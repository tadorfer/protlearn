import pytest
import numpy as np
from ..lasso import lasso

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_lasso():
    "Test Lasso-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute Lasso-based feature importances
    X_reduced = lasso(X, y)

    assert X.shape > X_reduced.shape