import pytest
import numpy as np
from ..moreau_broto import moreau_broto
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_moreau_broto():
    "Test Moreau-Broto Autocorrelation"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # test Moran's I
    mb = moreau_broto(X_list)

    # test array contents
    np.testing.assert_almost_equal(mb[0], np.array([
         0.30741363,  0.14572374,  0.24889737,  0.17774769,  0.31419645,
         0.59789614,  0.52090645,  0.54431542]))
    np.testing.assert_almost_equal(mb[1], np.array([
        -0.22151544,  0.48430045,  0.40262559,  0.65978741,  0.58163617,
         0.48620151,  1.57705014, -0.13977796]))
    np.testing.assert_almost_equal(mb[2], np.array([
         0.04011865, -0.05581411,  0.54491746,  0.48283298,  0.40420316,
         0.7053263 ,  0.07440376, -0.1527529]))

    # test ValueError (alphabetical)
    with pytest.raises(ValueError):
        mb_err = moreau_broto(X_err)

    # test ValueError (maximum lag)
    with pytest.raises(ValueError):
        mb_err = moreau_broto(X_list, d= 31) 

    # test ValueError (d >= min_len)
    with pytest.raises(ValueError):
        mb_err = moreau_broto(X_err, d=7)