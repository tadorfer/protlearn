import pytest
import numpy as np
from ..mutual_information import mutual_information

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_mutual_information():
    "Test MI-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute MI-based feature importances
    X_reduced = mutual_information(X, y, top=10)

    # test array shape
    X_reduced.shape == (700, 10)