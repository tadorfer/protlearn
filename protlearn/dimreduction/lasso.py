# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression

def lasso(X, y, C=1.0):
    """Select top features based on Lasso regularization.
    
    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
    
    y : ndarray of shape (n_samples,)
    
    C : float, default=1.0
        Inverse of regularization strength.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_features_post)
    
    """  
    
    mdl = SelectFromModel(LogisticRegression(C=C, 
                                             penalty='l1', 
                                             solver='liblinear'))
    arr = mdl.fit_transform(X, y)
        
    return arr