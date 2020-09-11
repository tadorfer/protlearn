# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from mlxtend.feature_selection import SequentialFeatureSelector

def sequential(X, y, estimator, direction='forward', n_features=10, cv=0):
    """Sequential feature selection.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 

    y : labels, ndarray of shape (n_samples,)
    
    estimator : object
        Classifier must include coef_ or feature_importances_ attribute.
        
    direction : string, default='forward'
        Direction of sequential model, can be 'forward' or 'backward'.
    
    n_features : int, default=None
        Number of features to select.
        
    cv : int, default=0
        Number of cross-validation steps.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features)
    
    """

    if direction == 'forward':
        method = True
    elif direction == 'backward':
        method = False
    
    mdl = SequentialFeatureSelector(estimator,
                                    k_features=n_features,
                                    forward=method,
                                    floating=False,
                                    verbose=0,
                                    scoring='accuracy',
                                    cv=cv)

    arr = mdl.fit_transform(X, y)
    
    return arr