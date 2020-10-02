# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier

def tree_importance(X, y, *, clf=None, method='random_forest', top=None, 
                    n_estimators=100, max_depth=None, importance_type='gain'):
    """Decision tree feature importance.

    This function returns the features that were selected as important by 
    decision tree algorithms such as Random Forest and XGBoost. 
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        Feature matrix.
    
    y : ndarray of shape (n_samples,)
        Response variables.

    clf : object or None, default=None
        Customized classifier.
    
    method : string, default='random_forest'
        'random_forest' : Random Forest Classifier
        'xgboost' : XGBoost Classifier
    
    top : int or None, default=None
        Number of top features to select.
        
    n_estimators : int or None, default=2
        Number of trees in the forest.
        
    max_depth : int or None, default=None
        Maximum depth of the tree.
        
    importance_type : string, default='gain'
        For XGBoost only:
        'gain' : average gain of splits which use the feature
        'weight' : number of times the a feature appears in the tree
        'cover' : average coverage of splits which use the feature
        'total_gain' : Total gain
        'total_cover' : Total cover

    Returns
    -------

    arr :  ndarray of shape (n_samples, top)
        Array containing the top features based on tree-importance.

    Examples
    --------

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import tree_importance
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = tree_importance(features, labels, top=10)
    >>> reduced.shape
    (3, 10)
    
    """
    
    if clf != None:
        clf = clf
    else:
        if method == 'random_forest':
            clf = RandomForestClassifier(n_estimators=n_estimators, 
                                        max_depth=max_depth)
        elif method == 'xgboost':
            clf = XGBClassifier(n_estimators=n_estimators, 
                                max_depth=max_depth, 
                                importance_type=importance_type)
        
    clf.fit(X, y)
    importances = clf.feature_importances_
    indices = np.argsort(-importances)[:top]
    arr = X[:,indices]
    
    return arr