import pytest
import numpy as np
from ..rf_importance import rf_importance

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_rf_importance():
    "Test RF-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute RF-based feature importances
    X_reduced = rf_importance(X, y, top=5)

    # test array shape
    X_reduced.shape == (700, 5)