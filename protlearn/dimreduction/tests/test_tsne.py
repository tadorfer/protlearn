import pytest
import numpy as np
from ..tsne import tsne

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_tsne():
    "Test t-SNE-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    X = X[:,:50]

    # perform t-SNE
    X_reduced = tsne(X, n_components=3, perplexity=5, pca_components=50)

    # test array shape
    X_reduced.shape == (700, 3)