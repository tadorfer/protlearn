import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protclass')
import numpy as np

from preprocessing import txt_to_df
from feature_engineering import length


def test_lengths():
    "Test sequence lengths"
    
    # load data
    df = txt_to_df(path+'/tests/docs/test_seq.txt', 0)
    
    # test integer lengths
    len_int = length(df, 'int')
    assert np.array_equal(len_int, np.array([6, 9, 7, 6]))
    
    # test one-hot-encoded lengths
    len_ohe = length(df, 'ohe')
    # columns: [6, 7, 9]
    assert np.array_equal(len_ohe, np.array([[1., 0., 0.],
                                             [0., 0., 1.],
                                             [0., 1., 0.],
                                             [1., 0., 0.]]))