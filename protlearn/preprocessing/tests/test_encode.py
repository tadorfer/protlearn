import pytest
import numpy as np
from preprocessing import encode
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_encode():
    "Test integer encoding"
    
   # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test encode
    enc_list, aa = encode(X_list)

    with pytest.raises(ValueError):
        enc_err, aa = encode(X_err)
    
    # test array contents
    assert np.array_equal(enc_list[0], np.array([1,  1, 15,  9, 20, 10, 10]))
    assert np.array_equal(enc_list[1], np.array([10,  4, 10,  2,  3, 13,  6, 13,  6]))
    assert np.array_equal(enc_list[2], np.array([15,  1,  1,  1, 12,  2,  3,  3]))

    # test padding
    enc_list, aa = encode(X_list, padding=True)
    assert enc_list.shape == (3, 9)
    assert [enc_list[0,i] == 0 for i in [7, 8]]
    assert [enc_list[2,i] == 0 for i in [8]]