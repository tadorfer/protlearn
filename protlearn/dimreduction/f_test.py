# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import SelectKBest, f_classif

def f_test(X, y, top=10):
    """Select top features based on ANOVA F-test scores.
    
    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
    
    y : ndarray of shape (n_samples,)
    
    top : int, default=10
        Choose number of top features to select.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_features_post)
    
    """  
    
    arr = SelectKBest(f_classif, k=top).fit_transform(X, y)
        
    return arr