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
    arr3, ng3 = ngram(X_list, n=3)

    # test ValueError (string input)
    with pytest.raises(ValueError):
        ngram_err, aa = ngram(X_err)

    # test ValueError (n > 3)
    with pytest.raises(ValueError):
        ngram_err, aa = ngram(X_err, n=4)

    # test array contents (n=2)
    np.testing.assert_almost_equal(arr2[0], np.array([
        0.16666667, 0.16666667, 0.16666667, 0.16666667, 0.16666667,
        0.16666667, 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.]))
    np.testing.assert_almost_equal(arr2[1], np.array([
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.125     , 0.125     , 0.125     , 0.125     ,
        0.125     , 0.125     , 0.25      , 0.        , 0.        ,
        0.        , 0.]))
    np.testing.assert_almost_equal(arr2[2], np.array([
        0.28571429, 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.14285714, 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.14285714, 0.14285714,
        0.14285714, 0.14285714]))

    # test array contents (n=3)
    np.testing.assert_almost_equal(arr3[0], np.array([
        0.2       , 0.2       , 0.2       , 0.2       , 0.2       ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.]))
    np.testing.assert_almost_equal(arr3[1], np.array([
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.14285714, 0.14285714, 0.14285714, 0.14285714, 0.14285714,
        0.14285714, 0.14285714, 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.]))
    np.testing.assert_almost_equal(arr3[2], np.array([
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.16666667, 0.16666667, 0.16666667,
        0.16666667, 0.16666667, 0.16666667]))

    # test absolute
    arr_rel, ng = ngram(X_list, n=2, method='absolute')
    assert np.array_equal(arr_rel[0], np.array([
        1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]))
    assert np.array_equal(arr_rel[1], np.array([
        0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 2., 0., 0., 0., 0.]))
    assert np.array_equal(arr_rel[2], np.array([
        2., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1.]))