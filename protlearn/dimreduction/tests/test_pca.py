import pytest
import numpy as np
from ..pca import pca

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_pca():
    "Test PCA-based dimensionality reduction"
    
    # load data
    features_w = np.load(PATH+'features.npy')
    features_large = np.load(PATH+'features_largeN.npy')

    # test warning (n_samples < n_features)
    features_reduced_w = pca(features_w)
    features_reduced_large = pca(features_large)

    # test array contents
    np.testing.assert_almost_equal(features_reduced_large[0,:], np.array([
         68.99258546, -81.64870178,  52.14384484,  42.11308778]), decimal=3)

    np.testing.assert_almost_equal(features_reduced_large[300,:], np.array([
         149.04942476, -95.54527396,  -6.43144111,  22.16476921]), decimal=3)

    np.testing.assert_almost_equal(features_reduced_large[-1,:], np.array([
         -5.88156252,  66.58515889, -48.50120568,  34.66457294]), decimal=3)

    # test array shape
    features_reduced_large.shape == (700, 10)