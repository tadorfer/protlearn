import pytest
import numpy as np
from ..moran import moran
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_moran():
    "Test Moran's I"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test Moran's I
    m = moran(X_list)

    # test array contents
    np.testing.assert_almost_equal(m[0], np.array([
         0.58182211,  0.23444196,  0.22908942,  0.14174852,  0.25901741,
         0.47640498,  0.60218565,  0.70541818]))
    np.testing.assert_almost_equal(m[1], np.array([
        -0.38994736,  0.26133405,  0.00794481,  0.41805449,  0.16466036,
         0.12089533,  0.43411552, -0.47697942]))
    np.testing.assert_almost_equal(m[2], np.array([
        -0.29862195, -0.05361181, -0.10350556,  0.41292954, -0.14182024,
        -0.11961891,  0.25280944, -0.49592705]))

    # test ValueError (alphabetical)
    with pytest.raises(ValueError):
        m_err = moran(X_err)

    # test ValueError (maximum lag)
    with pytest.raises(ValueError):
        m_err = moran(X_list, d= 31) 

    # test ValueError (d >= min_len)
    with pytest.raises(ValueError):
        m_err = moran(X_err, d=7)