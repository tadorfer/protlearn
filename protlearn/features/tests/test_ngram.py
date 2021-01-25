import pytest
import numpy as np
from ..ngram import ngram
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_ngram():
    "Test ngram composition"

    # load data
    X_list = open(PATH+'multiple.txt').read().splitlines()
    X_err = 'AGT2HT9'

    # test ngram (n=2, n=3)
    arr2, ng2 = ngram(X_list, n=2)
    assert arr2.shape == (3, 400)
    assert len(ng2) == 400
    arr3, ng3 = ngram(X_list, n=3)
    assert arr3.shape == (3, 8000)
    assert len(ng3) == 8000

    # test ValueError (string input)
    with pytest.raises(ValueError):
        ngram_err, aa = ngram(X_err)

    # test ValueError (n > 3)
    with pytest.raises(ValueError):
        ngram_err, aa = ngram(X_err, n=4)

    
    
    # test absolute
    arr_rel, ng = ngram(X_list, n=2, method='absolute')
    assert arr_rel.shape == (3, 400)
    assert len(ng) == 400
    