# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif, chi2, mutual_info_classif
from sklearn.preprocessing import MinMaxScaler

def univariate_filter(X, y, method='f_test', top=10):
    """Select top features based on univariate filtering.
    
    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
    
    y : ndarray of shape (n_samples,)
    
    method : string, default='f_test'
        'f_test' : ANOVA f-scores
        'chi2' : Chi-squared statistics
        'mutual_info' : Mutual information
    
    top : int, default=10
        Choose number of top features to select.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, top)
    
    """  
    
    if method == 'f_test':
        arr = SelectKBest(f_classif, k=top).fit_transform(X, y)
    elif method == 'chi2':
        # only non-negative features
        try:
            arr = SelectKBest(chi2, k=top).fit_transform(X, y)
        except:
            X = MinMaxScaler().fit_transform(X)
            arr = SelectKBest(chi2, k=top).fit_transform(X, y)
    elif method == 'mutual_info':
        arr = SelectKBest(mutual_info_classif, k=top).fit_transform(X, y)
    
    return arr