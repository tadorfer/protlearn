import numpy as np

from ..preprocessing import integer_encode


def test_integer_encode():
    "Test integer encoding"
    
    # load data
    data = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD']
    enc = integer_encode(data)
    
    # test array shape and type
    assert enc.shape == (4,)
    assert type(enc) == np.ndarray
    
    # test array contents
    assert np.array_equal(enc[0], np.array([1, 6, 17, 20, 10, 9]))
    assert np.array_equal(enc[1], np.array([18, 2, 8, 11, 11, 11, 13, 5, 13]))
    assert np.array_equal(enc[2], np.array([10, 15, 16, 1, 7, 7, 12]))
    assert np.array_equal(enc[3], np.array([1, 14, 4, 4, 19, 3]))

    # test padding
    enc = integer_encode(data, padding=True)
    assert enc.shape == (4, 9)
    assert [enc[0][i] == 0 for i in [6, 7, 8]]
    assert [enc[2][i] == 0 for i in [7, 8]]
    assert [enc[3][i] == 0 for i in [6, 7, 8]]