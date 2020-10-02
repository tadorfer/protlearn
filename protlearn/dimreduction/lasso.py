# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression

def lasso(X, y, C=1.0):
    """Lasso (L1) regularization.
    
    Linear Model trained with L1 prior as regularizer. 

    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
        Feature matrix.
    
    y : ndarray of shape (n_samples,)
        Response variables.
    
    C : float, default=1.0
        Inverse of regularization strength.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_features_post)
        Array containing lasso-reduced features.

    Examples
    --------

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import lasso
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = lasso(features, labels)
    >>> reduced.shape
    (3, 2)
    
    """  
    
    mdl = SelectFromModel(LogisticRegression(C=C, 
                                             penalty='l1', 
                                             solver='liblinear'))
    arr = mdl.fit_transform(X, y)
        
    return arr