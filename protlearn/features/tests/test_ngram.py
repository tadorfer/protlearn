import pytest
import numpy as np
from features import ngram
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_ngram():
    "Test ngram composition"

    # load data
    data_list = ['AGTY', 'AHHN', 'AQEE']
    data_str = 'AGTYLK'
    data_fasta = PATH+'sarcolipin.fasta'
    data_error = 'AGT2HT9'
    arr, ng = ngram(data_list, method='absolute')
    arr_str, ng = ngram(data_str, method='absolute')
    arr_fasta, ng = ngram(data_fasta, method='absolute')

    with pytest.raises(ValueError):
        enc_error, aa = ngram(data_error)

    # test list data
    assert np.array_equal(arr[0], np.array([1.,1.,1.,0.,0.,0.,0.,0.,0.]))
    assert np.array_equal(arr[1], np.array([0.,0.,0.,1.,1.,1.,0.,0.,0.]))
    assert np.array_equal(arr[2], np.array([0.,0.,0.,0.,0.,0.,1.,1.,1.]))

    # test string data
    assert np.array_equal(arr_str, np.ones((1,5)))

    # test fasta data 
    assert np.array_equal(arr_fasta, np.ones((1,30)))

