# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import RFE

def rfe(X, y, *, estimator, n_features=None, step=1):
    """Recursive feature elimination.

    This function selects features by recursively considering smaller and 
    smaller feature subsets. First, the estimator is trained on the initial 
    feature matrix and the importance of each feature is obtaibed through a 
    coef_ or a feature_importances_ attribute. Subsequently, the least 
    important features are pruned from the current feature subset. This is 
    repeated recursively on the pruned subset until the desired number of 
    features is eventually reached.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        Feature matrix.

    y : labels, ndarray of shape (n_samples,)
        Response variables.
    
    estimator : object
        Classifier - must include coef_ or feature_importances_ attribute.
        
    n_features : int or None, default=None
        Number of features to select. If None, half of the features are selected.
        
    step : int, default=1
        Number of features to remove at each iteration.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features)
        Array containing the RFE-selected features.
    
    ranking : ndarray of shape (n_features_pre,)
        Ranking of the features (with 1 being the best).

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import rfe
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> rf = RandomForestClassifier()
    >>> reduced, _ = rfe(features, labels, rf, n_features=10, step=5)
    >>> reduced.shape
    (3, 10)
    
    """
    
    selector = RFE(estimator, n_features_to_select=n_features, step=step)
    selector = selector.fit(X, y)
    arr = selector.transform(X)
    
    return arr, selector.ranking_