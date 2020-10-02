# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def tsne(X, *, n_components=2, perplexity=30, prior_pca=True, pca_components=50):
    """t-distributed stochastic neighbor embedding.

    t-SNE converts similarities between data points to joint probabilities and 
    tries to minimize the Kullback-Leibler divergence between the joint 
    probabilities of the low-dimensional embedding and the high-dimensional data.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        Feature matrix.
        
    n_components : int or None, default=2
        Dimension of embedded space.
        
    perplexity : int, default=30
        Related to the number of nearest neighbors that is used in other 
        manifold learning algorithms. Should be between 5 and 50. Larger 
        datasets require larger perplexity.
        
    prior_pca : bool, default=True
        It is recommended to reduce dimensionality before running t-SNE to 
        decrease computation time and noise.
        
    pca_components : int, default=50
        Dimension of PCA-preprocessed data that will serve as input to t-SNE.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_components)
        Array containing the t-SNE-transformed features.

    Examples
    --------

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import tsne
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = tsne(features, pca_components=3)
    >>> reduced.shape
    (3, 2)
    
    """
    
    if prior_pca:
        pca = PCA(n_components=pca_components)
        pca_reduced = pca.fit_transform(X)
        arr = TSNE(n_components=n_components, 
                   perplexity=perplexity).fit_transform(pca_reduced)
    else:
        arr = TSNE(n_components=n_components, 
                   perplexity=perplexity).fit_transform(X)
    
    return arr