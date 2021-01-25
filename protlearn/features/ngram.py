# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from itertools import product
from ..utils.validation import check_input, check_alpha, check_natural

def ngram(X, *, n=2, method='relative', start=1, end=None):
    """N-gram composition.
    
    This function computes the di- or tripeptide composition of amino acid 
    sequences. Therefore, the function argument 'n' can only take on 
    the values 2 and 3 - otherwise, it will raise a ValueError.
    
    Parameters
    ----------
    
    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
       
    n : int, default=2
        Integer denoting the desired n-gram composition.
        2 : dipeptide composition
        3 : tripepitde composition
        
    method : string, default='relative'
        'absolute': absolute n-gram composition
        'relative': relative n-gram composition

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_combinations)
        Depending on n, the returned array will be of size:
        - (n_samples, 400) for dipeptide composition
        - (n_samples, 8000) for tripeptide composition

    n-grams : list of length 400 or 8000
        List of n-grams corresponding to columns in arr.

    Examples
    --------

    >>> from protlearn.features import ngram
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> di, ngrams = ngram(seqs, n=2)
    >>> di.shape
    (2, 400)
    >>> len(ngrams)
    400
    >>> tri, ngrams = ngram(seqs, n=3)
    >>> tri.shape
    (2, 8000)
    >>> len(ngrams)
    8000

    """
    
    # input handling
    X = check_input(X)
    
    # make sure ngram is between 2-3
    valid = [2, 3]
    if n not in valid:
        raise ValueError("n must be one of %r." % valid)
        
    # compute n-gram combinations
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    combo = [''.join(i) for i in product(amino_acids, repeat=n)]
    combo_dict = {combo[i]: i for i in range(len(combo))}
    
    # initialize empty array
    arr = np.zeros((len(X), 20**n))
    
    # compute n-gram composition
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical
        check_natural(seq) # check for unnatural amino acids  
        seq = seq[start-1:end] # positional information
        for j in range(len(seq)-(n-1)):
            tmp = combo_dict[seq[j:j+n]]
            arr[i, tmp] += 1
                    
    if method=='absolute':
        return arr, combo

    elif method=='relative':
        for i in range(arr.shape[0]):
            arr[i,:] = arr[i,:]/sum(arr[i,:])
        return arr, combo