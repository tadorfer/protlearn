# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import SelectKBest, mutual_info_classif

def mutual_information(X, y, n_neighbors=3, top=10):
    """Select top features based on mutual information scores.
    
    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
    
    y : ndarray of shape (n_samples,)
    
    n_neighbors: int, default=3
        Number of neighbors to use for MI estimation for continuous variables.
    
    top : int, default=10
        Choose number of top features to select.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_features_post)
    
    """  
    
    arr = SelectKBest(mutual_info_classif, k=top).fit_transform(X, y)
        
    return arr