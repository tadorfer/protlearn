# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def lda(X, y, *, solver='svd', shrinkage=None, n_components=None):
    """Linear discriminant analysis.

    This function reduces the dimensionality of the input by projecting it to 
    the most discriminative directions.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        Feature matrix. 
    
    y : ndarray of shape (n_samples,)
        Response variables.
    
    solver : string, default='svd'
        'svd' : Singular value decomposition
        'lsqr' : Least squares solution
        'eigen' : Eigenvalue decomposition
        
    shrinkage : string, float, or None, default=None
        Shrinkage parameter.
        None : no shrinkage
        'auto' : automatic shrinkage using the Ledoit-Wolf lemma
        float between 0 and 1: fixed shrinkage parameter
        
    n_components : int or None, default=None
        Number of components for dimensionality reduction. This parameter 
        cannot be larger than min(n_features, n_classes - 1).

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features_post)
        Array containing the LDA-transformed features.

    Examples
    --------

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import lda
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = lda(features, labels, n_components=1)
    >>> reduced.shape
    (3, 1)
    
    """
    
    mdl = LinearDiscriminantAnalysis(solver=solver, 
                                     shrinkage=shrinkage,
                                     n_components=n_components)
    arr = mdl.fit_transform(X, y)
    
    return arr