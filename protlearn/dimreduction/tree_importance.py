# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier

def tree_importance(X, y, method='random_forest', top=None, n_estimators=100, 
                    max_depth=None, importance_type='gain'):
    """Dimensionality reduction using tree-based feature importances.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
    
    y : ndarray of shape (n_samples,)
    
    method : string, default='random_forest'
        'random_forest' : Random Forest Classifier
        'xgboost' : XGBoost Classifier
    
    top : int or None, default=None
        Top number of features to select.
        
    n_estimators : int or None, default=2
        Number of trees in the forest.
        
    max_depth : int or None, default=None
        Maximum depth of the tree.
        
    mportance_type : string, default='gain'
        For XGBoost only.
    
        'gain' : average gain of splits which use the feature
        'weight' : number of times the a feature appears in the tree
        'cover' : average coverage of splits which use the feature
        'total_gain' : Total gain
        'total_cover' : Total cover

    Returns
    -------

    arr :  ndarray of shape (n_samples, top)
    
    """
    
    if method == 'random_forest':
        mdl = RandomForestClassifier(n_estimators=n_estimators, 
                                     max_depth=max_depth)
    elif method == 'xgboost':
        mdl = XGBClassifier(n_estimators=n_estimators, 
                            max_depth=max_depth, 
                            importance_type=importance_type)
        
    mdl.fit(X, y)
    importances = mdl.feature_importances_
    indices = np.argsort(-importances)[:top]
    arr = X[:,indices]
    
    return arr