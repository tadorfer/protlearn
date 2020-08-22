import pytest
import numpy as np
from ..validation import check_input
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_validation():
    "Test validation function"

    # load data
    X_str = 'AGTYLK'
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_fasta_single = PATH+'sarcolipin.fasta'
    X_fasta_multiple = PATH+'multiple.fasta'
    X_err = 3975
    
    # test check_input
    val_str = check_input(X_str)
    val_list = check_input(X_list)
    val_fasta_single = check_input(X_fasta_single)
    val_fasta_multiple = check_input(X_fasta_multiple)

    # test TypeError
    with pytest.raises(TypeError):
        val_err = check_input(X_err)