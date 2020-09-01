import pytest
import numpy as np
from ..correlation import correlation

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_correlation():
    "Test correlation-based dimensionality reduction"
    
    # load data
    features = np.load(PATH+'features.npy')

    # test correlation
    corr = correlation(features)

    # test array contents
    np.testing.assert_almost_equal(corr, np.array([
        [0.28571429, 0.        , 0.        , 0.79428571],
        [0.        , 0.11111111, 0.11111111, 0.45      ],
        [0.375     , 0.125     , 0.25      , 0.60125   ]]), decimal=3)