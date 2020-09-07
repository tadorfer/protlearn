# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from xgboost import XGBClassifier

def xgb_importance(X, y, top=None, max_depth=6, importance_type='gain'):
    """Dimensionality reduction using XGBoost-based feature importances.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
    
    y : ndarray of shape (n_samples,)
    
    top : int or None, default=None
        Top number of features to select.
        
    max_depth : int, default=6
        Maximum depth of the tree.
        
    importance_type : string, default='gain'
    
        'gain' : average gain of splits which use the feature
        'weight' : number of times the a feature appears in the tree
        'cover' : average coverage of splits which use the feature
        'total_gain' : Total gain
        'total_cover' : Total cover

    Returns
    -------

    arr :  ndarray of shape (n_samples, top)
    
    """
    
    xgb = XGBClassifier(max_depth=max_depth, importance_type=importance_type)
    xgb.fit(X, y)
    importances = xgb.feature_importances_
    indices = np.argsort(-importances)[:top]
    arr = X[:,indices]
    
    return arr