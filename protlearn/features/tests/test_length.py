import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protlearn')
import numpy as np

from ..features import length


def test_lengths():
    "Test sequence lengths"
    
    # load data
    data = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD']
    
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