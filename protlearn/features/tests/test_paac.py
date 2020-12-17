import pytest
import numpy as np
from ..paac import paac
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_paac():
    "Test pseudo amino acid composition"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get paac
    paac_list, aa = paac(X_list, lambda_=1, remove_zero_cols=True)
    
    # test paac
    np.testing.assert_almost_equal(paac_list, 
    np.array([[1.79965836, 0.        , 0.        , 0.        , 0.        ,
        0.89982918, 1.79965836, 0.        , 0.        , 0.89982918,
        0.89982918, 0.10017082],
       [0.        , 0.92749507, 0.92749507, 0.92749507, 1.85499014,
        0.        , 1.85499014, 0.        , 1.85499014, 0.        ,
        0.        , 0.07250493],
       [2.77372003, 0.92457334, 1.84914669, 0.        , 0.        ,
        0.        , 0.        , 0.92457334, 0.        , 0.92457334,
        0.        , 0.07542666]])
    , decimal=3)

    # test ValueError
    with pytest.raises(ValueError):
        paac_error, aa = paac(X_err)