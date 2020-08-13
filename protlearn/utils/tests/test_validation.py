import pytest
import numpy as np
from utils import validation
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_validation():
    "Test validation function"

    # load data
    data_list = ['AGTY', 'AHHN', 'AQEE']
    data_str = 'AGTYLK'
    data_fasta_single = PATH+'sarcolipin.fasta'
    data_fasta_multiple = PATH+'multiple.fasta'
    data_error = 'AGT2HT9'
    
    # test check_input
    X_list = check_input(data_list)
    X_str = check_input(data_str)
    X_fasta_single = check_input(data_fasta_single)
    X_fasta_multiple = check_input(data_fasta_multiple)

    # test ValueError
    with pytest.raises(ValueError):
        X_err = check_input(data_error)