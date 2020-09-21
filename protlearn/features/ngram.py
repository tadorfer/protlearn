# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
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
    
    arr : ndarray of shape (n_samples, n_unique^n)
        Depending on n, the returned array will be of size:
        - (n_samples, 400) for dipeptide composition
        - (n_samples, 8000) for tripeptide composition
        if all possible n-gram combinations are represented.

    n-grams : list of length n_unique^n
        List of n-grams corresponding to columns in arr.

    Examples
    --------

    >>> from protlearn.features import ngram
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> di, ngrams = ngram(seqs, n=2)
    >>> di
    array([[0.25      , 0.25      , 0.25      , 0.25      , 0.        ,
            0.        , 0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.16666667, 0.16666667,
            0.16666667, 0.16666667, 0.16666667, 0.16666667]])
    >>> ngrams
    ['AR', 'KL', 'LY', 'RK', 'EE', 'ER', 'GL', 'KP', 'PG']
    >>> tri, ngrams = ngram(seqs, n=3)
    array([[0.33333333, 0.33333333, 0.33333333, 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.2       , 0.2       ,
            0.2       , 0.2       , 0.2       ]])
    >>> ngrams
    ['ARK', 'KLY', 'RKL', 'EER', 'ERK', 'KPG', 'PGL', 'RKP']

    """
    
    # input handling
    X = check_input(X)
    
    # make sure ngram is between 2-3
    valid = [2, 3]
    if n not in valid:
        raise ValueError("n must be one of %r." % valid)
    
    # compute n-gram composition
    ngdict = dict()
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical
        check_natural(seq) # check for unnatural amino acids   
        seq = seq[start-1:end]
        keys = [seq[i:i+n] for i in range(len(seq)-n+1)]
        unq = sorted(set(keys))
        ngram = sorted(Counter(keys).items())
        vals = [i[1] for i in ngram]
        for num, j in enumerate(unq):
            if j in ngdict:
                ngdict[j].append(vals[num])
            else:
                # put zeros before if present for first time
                ngdict[j] = i*[0]+[vals[num]]

        # append values not present in ngram with zero
        if i != 0:
            maxlen = max([len(c) for c in ngdict.values()])
            for z in ngdict.values():
                if len(z) < maxlen:
                    z.append(0)

    arr = np.array(list(ngdict.values()), dtype=float).T
    ngrams = list(ngdict.keys())
                    
    if method=='absolute':
        return arr, ngrams

    elif method=='relative':
        for i in range(arr.shape[0]):
            arr[i,:] = arr[i,:]/sum(arr[i,:])
        return arr, ngrams