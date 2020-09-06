import pytest
import numpy as np
from ..tsne import tsne

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_tsne():
    "Test t-SNE-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')

    # perform t-SNE
    X_reduced = tsne(X, n_components=3, perplexity=5, pca_components=50)

    # test array contents
    np.testing.assert_almost_equal(X_reduced[0,:], np.array([
         9.4649725, -6.909212 , -6.9694576]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[300,:], np.array([
         -20.854897  ,  -5.7361755 ,   0.21855178]), decimal=3)

    np.testing.assert_almost_equal(X_reduced[-1,:], np.array([
         3.7729113, -8.679973 , -7.2757826]), decimal=3)

    # test array shape
    X_reduced.shape == (700, 3)