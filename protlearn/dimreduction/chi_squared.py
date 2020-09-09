# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import SelectKBest, chi2

def chi_squared(X, y, top=10):
    """Select top features based on chi-squared scores.
    
    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
    
    y : ndarray of shape (n_samples,)
    
    top : int, default=10
        Choose number of top features to select.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_features_post)
    
    Notes
    -----
    
    If features contain negative values, minmax scaling will be performed.
    
    """  
    
    try:
        arr = SelectKBest(chi2, k=top).fit_transform(X, y)
    except:
        X = MinMaxScaler().fit_transform(X)
        arr = SelectKBest(chi2, k=top).fit_transform(X, y)
        
    return arr