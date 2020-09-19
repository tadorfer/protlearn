import pytest
import numpy as np
from ..aaindex1 import aaindex1
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_aaindex1():
    "Test AAIndex1"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get aaindex1
    aaind1, desc = aaindex1(X_list)
    
    # test shape
    assert aaind1.shape == (3, 553)
    
    # test some indices
    ANDN920101 = np.array([4.34, 4.31777778, 4.54375]) # index 0
    QIAN880126 = np.array([-.157, .019, .029]) # index 277
    KARS160122 = np.array([1.91385714, 4.55744444, 2.39225]) # index -1
    np.testing.assert_almost_equal(aaind1[:,0], ANDN920101, decimal=3)
    ind = np.where(desc=='QIAN880126')[0][0]
    np.testing.assert_almost_equal(aaind1[:,ind], QIAN880126, decimal=3)
    np.testing.assert_almost_equal(aaind1[:,-1], KARS160122, decimal=3)
    
    # test standardization (zscore)
    aaind1_z, desc = aaindex1(X_list, standardize='zscore')
    # test mean = 0
    for i in range(aaind1_z.shape[0]):
        assert abs(round(aaind1_z[:,1].mean())) == 0
    # test std --> 1
    for i in range(aaind1_z.shape[0]):
        assert round(aaind1_z[:,i].std(), 1) ==\
               round(aaind1_z[:,0].std(), 1)
        
    # test standardization (minmax)
    aaind1_mm, desc = aaindex1(X_list, standardize='minmax')
    # test minimum and maximum
    for i in range(aaind1_mm.shape[0]):
        assert round(aaind1_mm[:,i].min()) == 0
        assert round(aaind1_mm[:,i].max()) == 1

    # test ValueError
    with pytest.raises(ValueError):
        aaind1_error, desc = aaindex1(X_err)