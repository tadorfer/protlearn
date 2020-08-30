import pytest
import numpy as np
from ..motif import motif
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_motif():
    "Test custom motif"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get motif
    motif_1 = motif(X_list, pattern='AA[RN]{K}')
    motif_2 = motif(X_list, pattern='xxxCD{D}')
    
    # test motif
    np.testing.assert_almost_equal(motif_1, np.array([0., 0., 1.]), decimal=3)
    np.testing.assert_almost_equal(motif_2, np.array([0., 1., 0.]), decimal=3)

    # test ValueError
    with pytest.raises(ValueError):
        motif_error = motif(X_err, pattern='AA[RN]{K}')