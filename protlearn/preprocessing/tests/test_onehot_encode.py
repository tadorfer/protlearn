import pytest
import numpy as np
from ..onehot_encode import onehot_encode
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_onehot_encode():
    "Test one-hot encoding"
    
    # load data
    X_str = 'AYTLG'
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test encode
    enc_str = onehot_encode(X_str)
    enc_list = onehot_encode(X_list)

    with pytest.raises(ValueError):
        enc_err = onehot_encode(X_err)
    
    # test OHE contents
    assert np.array_equal(enc_str, np.array([[
        [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 1.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         1., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0.]]]))
    

    # test OHE shape
    assert enc_str.shape == (1, 5, 20)
    assert enc_list.shape == (3, 9, 20)