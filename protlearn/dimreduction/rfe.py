# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import RFE

def rfe(X, y, estimator, n_features=None, step=1):
    """Feature ranking with recursive feature elimination.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 

    y : labels, ndarray of shape (n_samples,)
    
    estimator : object
        Classifier must include coef_ or feature_importances_ attribute.
        
    n_features : int or None, default=None
        Number of features to select. If None, half of the features are selected.
        
    step : int, default=1
        Number of features to remove at each iteration.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features)
    
    ranking : ndarray of shape (n_features_pre,)
    
    """
    
    selector = RFE(estimator, n_features_to_select=n_features, step=step)
    selector = selector.fit(X, y)
    arr = selector.transform(X)
    
    return arr, selector.ranking_