# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def lda(X, y, solver='svd', shrinkage=None, n_components=None):
    """Dimensionality reduction using Linear Discriminant Analysis.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
    
    y : ndarray of shape (n_samples,)
    
    solver : string, default='svd'
        'svd' : Singular value decomposition
        'lsqr' : Least squares solution
        'eigen' : Eigenvalue decomposition
        
    shrinkage : string, float, or None, default=None
        None : no shrinkage
        'auto' : automatic shrinkage using the Ledoit-Wolf lemma
        float between 0 and 1: fixed shrinkage parameter
        
    n_components : int or None, default=None
        Number of components for dimensionality reduction.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_features_post)
    
    """
    
    mdl = LinearDiscriminantAnalysis(solver=solver, 
                                     shrinkage=shrinkage,
                                     n_components=n_components)
    arr = mdl.fit_transform(X, y)
    
    return arr