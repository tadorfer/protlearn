# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
from ..utils.validation import check_input, check_alpha, check_natural



def posrich(X, *, position, aminoacid):
    """Position-specific amino acids.

    This function returns a binary vector or matrix in which ones indicate the 
    presence of the given amino acid(s) at the specified position(s), and zeros
    indicate their absence. 

    Parameters
    ----------
    
    X : string, fasta, or a list thereof
        Dataset of amino acid sequences.
       
    position : int or list
        Integer or list of integers denoting the position(s) in the sequence. 

    aminoacid : string or list
        String or list of strings indicating the amino acid(s) of interest.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, ) or (n_samples, n_positions)
        Binary vector/matrix indicating position-specific presence of amino acids.

    Notes
    -----

    The position argument is based on one-based indexing.

    Examples
    --------

    >>> from protlearn.features import posrich
    >>> seqs = ['ARKLY', 'ERNLAPG', 'YRLQLLLY']   
    >>> pos_single = posrich(seqs, position=4, aminoacid='L')
    >>> pos_single
    array([1., 1., 0.])
    >>> pos_multiple = posrich(seqs, position=[2,3,4], aminoacid=['R','N','L'])
    array([[1., 0., 1.],
           [1., 1., 1.],
           [1., 0., 0.]])
       
    """
    
    # input handling
    X = check_input(X)
    
    if isinstance(position, int) and isinstance(aminoacid, str):
        arr = np.zeros((len(X),))
        for a, seq in enumerate(X):
            if str.isalpha(seq) == True:
                pass
            else:
                raise ValueError('Data type must be string!')
            for i, aa in enumerate(seq):
                if i == position-1 and aa == aminoacid:
                    arr[a] = 1
        return arr
    
    elif isinstance(position, list) and isinstance(aminoacid, list):
    
        if len(position) != len(aminoacid):
            raise ValueError("Number of positions does not match number of amino acids")

        arr = np.zeros((len(X), len(position)))
        for a, seq in enumerate(X):
            check_alpha(seq) # check if alphabetical  
            check_natural(seq) # check for unnatural amino acids
            for i in range(len(position)):
                if seq[position[i]-1] == aminoacid[i]:
                    arr[a, i] = 1
        return arr
    
    else:
        raise ValueError("The arguments position and aminoacid must either be integer/string or lists of integers/strings.")