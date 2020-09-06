# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def tsne(X, n_components=2, perplexity=30, prior_pca=True, pca_components=50):
    """Dimensionality reduction using t-distributed stochastic neighbor embedding.
    
    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features_pre) 
        
    n_components : int or None, default=2
        Dimension of embedded space.
        
    perplexity : int, default=30
        Related to the number of nearest neighbors that is used in other 
        manifold learning algorithms. Should be between 5 and 50. Larger 
        datasets require larger perplexity.
        
    prior_pca : bool, default=True
        It is recommended to reduce dimensionality before running t-SNE to 
        decrease computation time and noise.
        
    pca_components: int, default=50
        Dimension of PCA-preprocessed data that will serve as input to t-SNE.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_components)
    
    """
    
    if prior_pca:
        pca = PCA(n_components=pca_components)
        pca_reduced = pca.fit_transform(X)
        pca_reduced
        arr = TSNE(n_components=n_components, 
                   perplexity=perplexity).fit_transform(pca_reduced)
    else:
        arr = TSNE(n_components=n_components, 
                   perplexity=perplexity).fit_transform(X)
    
    return arr