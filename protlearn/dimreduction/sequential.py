# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from mlxtend.feature_selection import SequentialFeatureSelector

def sequential(X, y, *, estimator, direction='forward', n_features=10, cv=0):
    """Sequential feature selection.

    Sequential feature selection algorithms are a family of greedy search 
    algorithms that are used to reduce an initial d-dimensional feature space 
    to a k-dimensional feature subspace where k < d. These algorithms remove or 
    add one feature at a time based on the classifier performance until a 
    feature subset of the desired size k is reached.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        Feature matrix.

    y : labels, ndarray of shape (n_samples,)
        Response variables.
    
    estimator : object
        Classifier - must include coef_ or feature_importances_ attribute.
        
    direction : string, default='forward'
        Direction of sequential model, can be 'forward' or 'backward'.
    
    n_features : int, default=None
        Number of features to select.
        
    cv : int, default=0
        Number of cross-validation steps.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features)
        Array containing features selected by the sequential models.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import sequential
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> rf = RandomForestClassifier()
    >>> reduced = sequential(features, labels, rf, n_features=10)
    >>> reduced.shape
    (3, 10)
    
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