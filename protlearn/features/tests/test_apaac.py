import pytest
import numpy as np
from ..apaac import apaac
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_apaac():
    "Test amphiphilic pseudo amino acid composition"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get apaac
    apaac_list, aa = apaac(X_list, lambda_=1, remove_zero_cols=True)
    
    # test apaac
    np.testing.assert_almost_equal(apaac_list, 
    np.array([[ 1.90085811,  0.        ,  0.        ,  0.        ,  0.        ,
         0.95042906,  1.90085811,  0.        ,  0.        ,  0.95042906,
         0.95042906,  0.03019396,  0.01937699],
       [ 0.        ,  1.02931611,  1.02931611,  1.02931611,  2.05863221,
         0.        ,  2.05863221,  0.        ,  2.05863221,  0.        ,
         0.        , -0.00987507, -0.01944104],
       [ 2.98153907,  0.99384636,  1.98769271,  0.        ,  0.        ,
         0.        ,  0.        ,  0.99384636,  0.        ,  0.99384636,
         0.        , -0.00718136,  0.01333501]])
    , decimal=3)

    # test ValueError
    with pytest.raises(ValueError):
        apaac_error, aa = apaac(X_err)