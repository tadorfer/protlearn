# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from ..utils.validation import check_input, check_alpha, check_natural

def entropy(X, *, standardize='none', start=1, end=None):
    """Shannon entropy.

    This function computes the Shannon entropy for each sequence in the 
    dataset.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.

    standardize : string, default='none'
        'none' : unstandardized matrix will be returned
        'zscore' : matrix is standardized to have
                   a mean of 0 and standard deviation of 1.
        'minmax' : matrix is normalized to have a range of [0, 1].

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 1) if len(X) > 1, otherwise float
        Array containing Shannon entropy values for each sequence.

    Examples
    --------

    >>> from protlearn.features import entropy
    >>> seqs = ['ARKLY', 'EERKPGL', 'AAAAAALY']
    >>> ent = entropy(seqs)
    >>> ent
    array([[2.32192809], [2.52164064], [0.64020643]])

   """ 
    
    # input handling
    X = check_input(X)
    
    # initialize empty array with shape (n_samples, 1)
    arr = np.zeros((len(X), 1))
    
    # compute shannon entropy
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        cnt = Counter(seq)
        ent_aa = {k: (cnt[k]/len(seq))*np.log2(cnt[k]/len(seq)) for k in cnt.keys()}
        arr[i] = -sum(ent_aa.values())
        
    if len(arr) == 1:
        return arr[0][0]
    
    if standardize == 'none':
        return arr

    elif standardize == 'zscore':
        scaler = StandardScaler().fit(arr)
        return scaler.transform(arr)
    
    elif standardize == 'minmax':
        scaler = MinMaxScaler().fit(arr)
        return scaler.transform(arr)