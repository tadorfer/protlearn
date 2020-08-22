# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
from ..utils.validation import check_input

def ngram(X, n=2, method='relative', start=1, end=None):
    """Compute n-gram peptide composition.
    
    This function computes the di- or tripeptide composition of amino acid 
    sequences. Therefore, the function argument 'n' can only take on 
    the values 2 and 3 - otherwise, it will raise a ValueError.
    
    Parameters
    ----------
    
    X : string, fasta, or a list thereof 
       
    n : int, default=2
        Integer denoting the desired n-gram composition.
        
        2 : dipeptide composition
        3 : tripepitde composition
        
    method: string, default='relative'
    
        'absolute': compute absolute ngram composition
        'relative': compute relative ngram composition

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, n_unique_20^ngram)
        Depending on ngram, the returned dataframe will be of size:
        - (n_samples, 400) for dipeptide composition
        - (n_samples, 8000) for tripeptide composition
        if all possible ngram combinations are represented.
       
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
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')

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