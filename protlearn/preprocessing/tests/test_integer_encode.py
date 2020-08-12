import numpy as np
from preprocessing import integer_encode
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_integer_encode():
    "Test integer encoding"
    
    # load data
    data_list = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD']
    data_str = 'AGTYLK'
    data_fasta = PATH+'sarcolipin.fasta'
    enc, aa = integer_encode(data_list)
    enc_str, aa = integer_encode(data_str)
    enc_fasta, aa = integer_encode(data_fasta)
    
    
    # test array contents
    assert np.array_equal(enc[0], np.array([1, 6, 17, 20, 10, 9]))
    assert np.array_equal(enc[1], np.array([18, 2, 8, 11, 11, 11, 13, 5, 13]))
    assert np.array_equal(enc[2], np.array([10, 15, 16, 1, 7, 7, 12]))
    assert np.array_equal(enc[3], np.array([1, 14, 4, 4, 19, 3]))

    # test padding
    enc, aa = integer_encode(data_list, padding=True)
    assert enc.shape == (4, 9)
    assert [enc[0,i] == 0 for i in [6, 7, 8]]
    assert [enc[2,i] == 0 for i in [7, 8]]
    assert [enc[3,i] == 0 for i in [6, 7, 8]]