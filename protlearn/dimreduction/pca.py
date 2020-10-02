# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import warnings
import numpy as np
from sklearn.decomposition import PCA

def pca(X, *, thres=.9, whiten=False):
    """Principal component analysis.

    PCA is defined as an orthogonal linear transformation that transforms the 
    data to a new coordinate system such that the greatest variance by some 
    scalar projection of the data comes to lie on the first coordinate (called 
    the first principal component), the second greatest variance on the second 
    coordinate, and so on.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre)
        Feature matrix. 

    thres : float, default=.9
        Specify the desired explained variance.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features_post)
        Array containing the PCA components comprising the specified variance.

    Notes
    -----

    For the output to be meaningful, the number of samples should be larger than
    the number of features.

    Examples
    --------

    >>> from protlearn.dimreduction import pca
    >>> features.shape #from a larger dataset (not shown here)
    (1000, 575)
    >>> reduced = pca(features, thres=.9)
    (1000, 32)

    """

    # check input dimensionality
    if X.shape[0] < X.shape[1]:
        warnings.warn("The number of samples (%i) is less than the number of "
                      "features (%i). Therefore, the PCA output may not be "
                      "meaningful." % (X.shape[0], X.shape[1]))
    
    # fit and transform PCA
    pca = PCA(whiten=whiten).fit(X)
    var = pca.explained_variance_ratio_[0]
    comp = 1
    while var <= thres:
        var += pca.explained_variance_ratio_[comp]
        comp += 1
    
    arr = pca.transform(X)
    
    return arr[:,:comp]