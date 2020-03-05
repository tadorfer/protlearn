import os
import sys
sys.path.insert(0, os.path.dirname(os.getcwd())+'/ProtClass/protclass')
print(sys.path)
import numpy as np

from preprocessing import txt_to_df
from preprocessing import integer_encode


def test_integer_encode():
    "Test integer encoding"
    
    # load data
    df = txt_to_df('docs/test_seq.txt', 0)
    enc, labels = integer_encode(df)
    
    # test array shape and type
    assert enc.shape == (3,)
    assert type(enc) == np.ndarray
    
    # test array contents
    assert np.array_equal(enc[0], np.array([1, 6, 17, 20, 10, 9]))
    assert np.array_equal(enc[1], np.array([18, 17,  1, 20,  9, 10, 10, 10,  5]))
    assert np.array_equal(enc[2], np.array([ 1, 15,  1, 10, 17, 20,  8]))
    
    # test labels
    assert len(labels) == 3
    
    for i in range(len(labels)):
        assert labels[i] == 0

    # test padding
    enc, labels = integer_encode(df, padding=True)
    assert enc.shape == (3, 9)
    assert [enc[0][i] == 0 for i in [6, 7, 8]]
    assert [enc[2][i] == 0 for i in [7, 8]]