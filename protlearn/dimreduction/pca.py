# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import warnings
import numpy as np
from sklearn.decomposition import PCA

def pca(X, thres=.9, whiten=False):
    """Compute PCA components with specific variance.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 

    cutoff: float, default=.9
        Specify the desired variance.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features_post)

    """

    # check input dimensionality
    if X.shape[0] < X.shape[1]:
        warnings.warn("The number of samples (%i) is less than the number of "
                      "features (%i). Therefore, the PCA output may not be "
                      "meaningful." % (X.shape[0], X.shape[1]))
    
    # fit and transform PCA
    pca = PCA(whiten).fit(X)
    var = pca.explained_variance_ratio_[0]
    comp = 1
    while var <= thres:
        var += pca.explained_variance_ratio_[comp]
        comp += 1
    
    arr = pca.transform(X)
    
    return arr[:,:comp]