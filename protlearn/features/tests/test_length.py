import pytest
import numpy as np
from features import length

def test_lengths():
    "Test sequence lengths"
    
    # load data
    data = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD']
    data_str = 'AGTYLK'
    data_error = 'AGT2HT9'
    
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

    # test string data
    len_str = length(data_str, 'int')
    assert len_str == 6
    
    # test ValueError
    with pytest.raises(ValueError):
        len_error = length(data_error)