# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from sklearn.preprocessing import OneHotEncoder
from ..utils.validation import check_input, check_alpha, check_natural

def length(X, *, method='int'):
    """Sequence length in amino acids.
    
    The number of amino acids that a protein or peptide is comprised of will be
    counted and returned as either an integer (single sequence), array of 
    integers, or a one-hot-encoded array.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.

    method : string, default='int'
        'int' : interval data 
        'ohe' : one-hot encoded data

    Returns
    -------

    arr : int or ndarray of shape (n_samples, 1) or (n_samples, n_unique_lengths) 
        Array containing sequence lengths.

    Examples
    --------

    >>> from protlearn.features import length
    >>> seqs = ['ARKLY', 'EERKPGL', 'LLYPGP']
    >>> l_int = length(seqs)
    >>> l_int
    array([[5], [7], [6]])
    >>> l_ohe = length(seqs, method='ohe')
    array([[1., 0., 0.],
           [0., 0., 1.],
           [0., 1., 0.]])

    """
    
    # input handling
    X = check_input(X)
    
    # compute lengths
    arr = np.zeros((len(X), 1))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        arr[i] = len(seq)

    if method == 'int':
        # for single sequence return integer
        if len(arr) == 1:
            return arr[0]
        # for multiple sequences return array
        else:
            return arr
    elif method == 'ohe':
        encoder = OneHotEncoder(sparse=False)
        arr_ohe = encoder.fit_transform(arr)
        return arr_ohe