# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from sklearn.preprocessing import OneHotEncoder
from utils.validation import check_input

def length(X, method='int'):
    """Compute the length of proteins or peptides.
    
    The number of amino acids that a protein or peptide is comprised of will be
    counted and returned as either an integer (single sequence), array of integers,
    or a one-hot-encoded array.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    method : string, default='int'

        'int' : interval data 
        'ohe' : one-hot encoded data

    Returns
    -------

    l : ndarray of shape (n_samples, ) if method = 'int'
        ndarray of shape (n_samples, n_unique_lengths) if method = 'ohe'
        int if len(X) == 1

    """
    
    # input handling
    X = check_input(X)
    
    # compute lengths
    l = []
    for seq in X:
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical.')
        l.append(len(seq))

    if method == 'int':
        # for single sequence return integer
        if len(l) == 1:
            return l[0]
        # for multiple sequences return array
        else:
            return np.asarray(l)
    elif method == 'ohe':
        l = np.asarray(l).reshape((len(l), 1))
        encoder = OneHotEncoder(sparse=False)
        l_ohe = encoder.fit_transform(l)
        return l_ohe