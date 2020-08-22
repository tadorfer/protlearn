# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
from ..utils.validation import check_input



def posrich(X, position, aminoacid):
    """Compute the presence of an amino acid at a specific position.

    This function returns a binary feature vector or matrix in which ones 
    indicate the presence of the given amino acid(s) at the specified 
    position(s), and zeros indicate their absence. 

    Parameters
    ----------
    
    X : string, fasta, or a list thereof
       
    position : int or list
        Integer or list of integers denoting the position(s) in the sequence.

    aminoacid : string or list
        String or list of strings indicating the amino acid(s) of interest.
        
    Returns
    -------
    
    pos : ndarray of shape (n_samples, ) or (n_samples, n_positions)
       
    """
    
    # input handling
    X = check_input(X)
    
    if isinstance(position, int) and isinstance(aminoacid, str):
        pos = np.zeros((len(X),))
        for a, seq in enumerate(X):
            if str.isalpha(seq) == True:
                pass
            else:
                raise ValueError('Data type must be string!')
            for i, aa in enumerate(seq):
                if i == position-1 and aa == aminoacid:
                    pos[a] = 1
        return pos
    
    elif isinstance(position, list) and isinstance(aminoacid, list):
    
        if len(position) != len(aminoacid):
            raise ValueError("Number of positions does not match number of amino acids")

        pos = np.zeros((len(X), len(position)))
        for a, seq in enumerate(X):
            if str.isalpha(seq) == True:
                pass
            else:
                raise ValueError('Data type must be alphabetical!')
                
            for i in range(len(position)):
                if seq[position[i]-1] == aminoacid[i]:
                    pos[a, i] = 1
        return pos
    
    else:
        raise ValueError("The arguments position and aminoacid must either be integer/string or lists of integers/strings.")