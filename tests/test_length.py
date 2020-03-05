import os
import sys
sys.path.insert(0, os.path.dirname(os.getcwd())+'/protclass')
import numpy as np

from preprocessing import txt_to_df, integer_encode
from feature_engineering import length


def test_lengths():
    "Test sequence lengths"
    
    # load data
    df = txt_to_df('docs/test_seq.txt', 0)
    
    # test integer lengths
    len_int = length(df, 'int')
    assert np.array_equal(len_int, np.array([6, 9, 7]))
    
    # test one-hot-encoded lengths
    len_ohe = length(df, 'ohe')
    assert np.array_equal(len_ohe, np.array([[1., 0., 0.],
                                             [0., 0., 1.],
                                             [0., 1., 0.]]))