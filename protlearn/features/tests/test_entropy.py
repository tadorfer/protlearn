import pytest
import numpy as np
from ..entropy import entropy 
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_entropy():
    "Test Shannon entropy"
    
    # load data
    X_str = 'AYTLG'
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test entropy function
    np.testing.assert_almost_equal(2.3219, entropy(X_str), decimal=3)
    ent_list = entropy(X_list)
    np.testing.assert_array_almost_equal(ent_list, \
        np.array([[2.2359], [2.5032], [2.1556]]), decimal=3)

    # test standardization
    ent_zscore = entropy(X_list, standardize='zscore')
    np.testing.assert_array_almost_equal(ent_zscore, \
        np.array([[-.4195], [1.3793], [-.9598]]), decimal=3)
    ent_minmax = entropy(X_list, standardize='minmax')
    np.testing.assert_array_almost_equal(ent_minmax, \
        np.array([[.2309], [1.0], [0.]]), decimal=3)
    
    # test ValueError
    with pytest.raises(ValueError):
        ent_err = entropy(X_err)