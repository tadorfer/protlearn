import pytest
import numpy as np
from xgboost import XGBClassifier
from ..tree_importance import tree_importance

import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'test_data/')

def test_tree_importance():
    "Test tree-based dimensionality reduction"
    
    # load data
    X = np.load(PATH+'features_largeN.npy')
    y = np.load(PATH+'features_largeN_labels.npy')

    # compute tree-based feature importances
    X_rf, _ = tree_importance(X, y, top=5)
    X_xgb, _ = tree_importance(X, y, method='xgboost', top=5)

    # test array contents
    np.testing.assert_almost_equal(X_xgb[0,:], np.array([
        -0.01,  1.14444444,  4.73777778,  0.16666667, -0.67777778]), decimal=3)

    np.testing.assert_almost_equal(X_xgb[300,:], np.array([
        0.82777778,  1.03333333,  6.01888889,  0.07333333, -0.57888889]), decimal=3)

    np.testing.assert_almost_equal(X_xgb[-1,:], np.array([
        -0.52222222,  0.95555556,  6.02444444,  0.82555556, -0.97555556]), decimal=3)

    # test customized classifier
    clf = XGBClassifier()
    X_clf, _ = tree_importance(X, y, clf=clf, top=5)

    np.testing.assert_almost_equal(X_clf[0,:], np.array([
        -0.01,  1.14444444,  4.73777778,  0.16666667, -0.67777778]), decimal=3)

    np.testing.assert_almost_equal(X_clf[300,:], np.array([
        0.82777778,  1.03333333,  6.01888889,  0.07333333, -0.57888889]), decimal=3)

    np.testing.assert_almost_equal(X_clf[-1,:], np.array([
        -0.52222222,  0.95555556,  6.02444444,  0.82555556, -0.97555556]), decimal=3)

    # test array shape
    X_rf.shape == (700, 5)
    X_xgb.shape == (700,5)