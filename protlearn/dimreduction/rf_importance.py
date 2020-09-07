# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from sklearn.ensemble import RandomForestClassifier

def rf_importance(X, y, top=None, n_estimators=100, max_depth=None):
    """Dimensionality reduction using Random Forest-based feature importances.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
    
    y : ndarray of shape (n_samples,)
    
    top : int or None, default=None
        Top number of features to select.
        
    n_estimators : int or None, default=2
        Number of trees in the forest.
        
    max_depth : int or None, default=None
        Maximum depth of the tree.

    Returns
    -------

    arr :  ndarray of shape (n_samples, top)
    
    """
    
    rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth)
    rf.fit(X, y)
    importances = rf.feature_importances_
    indices = np.argsort(-importances)[:top]
    arr = X[:,indices]
    
    return arr