import pytest
import numpy as np
from ..correlation import correlation
from ...features.aac import aac
from ...features.aaindex1 import aaindex1
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_correlation():
    "Test correlation-based dimensionality reduction"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    
    # compute features
    aac = aac(X_list)[0]
    aaindex1 = aaindex1(X_list)[0]
    features = np.concatenate([aac, aaindex1], axis=1)

    # test correlation
    corr = correlation(X_list, cutoff=)

    # test array contents
    np.testing.assert_almost_equal(corr, np.array([
        [0.28571429, 0.        , 0.        , 0.79428571],
        [0.        , 0.11111111, 0.11111111, 0.45      ],
        [0.375     , 0.125     , 0.25      , 0.60125   ]])