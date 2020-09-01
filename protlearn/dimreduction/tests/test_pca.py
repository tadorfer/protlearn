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
         68.99258546, -81.64870178,  52.14384484,  42.11308778,
         0.363062  ,  21.56867107,  -3.07921867, -10.38377288,
         1.05282693, -13.17569411]), decimal=3)

    np.testing.assert_almost_equal(features_reduced_large[300,:], np.array([
         149.04942476, -95.54527396,  -6.43144111,  22.16476921,
         5.71360977,   6.63136857, -15.1653833 ,   7.76967884,
        -0.1647698 ,   7.37062579]), decimal=3)

    np.testing.assert_almost_equal(features_reduced_large[-1,:], np.array([
         -5.88156252,  66.58515889, -48.50120568,  34.66457294,
         -4.01704368,  -0.38098494,  -0.35958021,  -4.5810833 ,
         -13.30715876,   2.17495775]), decimal=3)

    # test array shape
    features_reduced_large.shape == (700, 10)