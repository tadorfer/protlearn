import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protlearn')
import numpy as np

from feature_engineering import length


def test_lengths():
    "Test sequence lengths"
    
    # load data
    data = open('/tests/docs/test_seq.txt', 'r').read().splitlines()
    
    # test integer lengths
    len_int = length(data, 'int')
    assert np.array_equal(len_int, np.array([6, 9, 7, 6]))
    
    # test one-hot-encoded lengths
    len_ohe = length(data, 'ohe')
    # columns: [6, 7, 9]
    assert np.array_equal(len_ohe, np.array([[1., 0., 0.],
                                             [0., 0., 1.],
                                             [0., 1., 0.],
                                             [1., 0., 0.]]))