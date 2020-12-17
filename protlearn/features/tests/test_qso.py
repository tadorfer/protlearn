import pytest
import numpy as np
from ..qso import qso
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_qso():
    "Test quasi-sequence-order"
    
    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'
    
    # get qso
    qso_sw, qso_g, desc = qso(X_list, d=1, remove_zero_cols=True)
    
    # test qso
    np.testing.assert_almost_equal(qso_sw, 
    np.array([[1.69679237, 0.        , 0.        , 0.        , 0.        ,
        0.84839618, 1.69679237, 0.        , 0.        , 0.84839618,
        0.84839618, 0.15160382],
       [0.        , 0.76052943, 0.76052943, 0.76052943, 1.52105887,
        0.        , 1.52105887, 0.        , 1.52105887, 0.        ,
        0.        , 0.23947057],
       [2.53298816, 0.84432939, 1.68865877, 0.        , 0.        ,
        0.        , 0.        , 0.84432939, 0.        , 0.84432939,
        0.        , 0.15567061]]),
        decimal=3)

    np.testing.assert_almost_equal(qso_g, 
    np.array([[9.19497954e-04, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
        0.00000000e+00, 4.59748977e-04, 9.19497954e-04, 0.00000000e+00,
        0.00000000e+00, 4.59748977e-04, 4.59748977e-04, 9.99540251e-01],
       [0.00000000e+00, 8.47644396e-05, 8.47644396e-05, 8.47644396e-05,
        1.69528879e-04, 0.00000000e+00, 1.69528879e-04, 0.00000000e+00,
        1.69528879e-04, 0.00000000e+00, 0.00000000e+00, 9.99915236e-01],
       [4.41748145e-04, 1.47249382e-04, 2.94498763e-04, 0.00000000e+00,
        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.47249382e-04,
        0.00000000e+00, 1.47249382e-04, 0.00000000e+00, 9.99852751e-01]]),
        decimal=3)

    # test ValueError (alphabetical)
    with pytest.raises(ValueError):
        qso_error_sw, qso_error_g, desc = qso(X_err, d=1)

    # test ValueError (exceeding maximum lag)
    with pytest.raises(ValueError):
        qso_error_sw, qso_error_g, desc = qso(X_err, d=15)