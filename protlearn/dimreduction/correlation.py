# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np

def correlation(X, thres=.9):
    """Pearson correlation.

    This function returns the features whose Pearson correlation with one 
    another is below a specified threshold, thus circumventing the problem of 
    multicollinearity.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        Feature matrix.

    thres : float, default=.9
        Features whose correlation coefficient is higher than this threshold 
        value will be removed.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features_post)
        Array containing features that correlate below the threshold with one 
        another.

    Examples
    --------

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import correlation
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = correlation(features, thres=.99)
    >>> reduced.shape
    (3, 12)

    """

    # remove correlated features
    corr = np.absolute(np.corrcoef(X, rowvar=False))
    upper = corr*np.triu(np.ones(corr.shape), k=1).astype(np.bool)
    to_drop = [column for column in range(upper.shape[1]) \
               if any(upper[:,column] > thres)]
    arr = np.delete(X, to_drop, axis=1)
    
    return arr