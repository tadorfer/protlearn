# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np

def correlation(X, cutoff=.9):
    """Remove highly correlated feature columns.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 

    cutoff: float, default=.9
        Remove features with correlation higher than cutoff value.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features_post)

    """

    # remove correlated features
    corr = np.absolute(np.corrcoef(X, rowvar=False))
    upper = corr*np.triu(np.ones(corr.shape), k=1).astype(np.bool)
    to_drop = [column for column in range(upper.shape[1]) \ 
               if any(upper[:,column] > cutoff)]
    arr = np.delete(X, to_drop, axis=1)
    
    return arr