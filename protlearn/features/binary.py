# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from ..utils.validation import check_input, check_alpha

def binary(X, padding=False, start=1, end=None):
    """Compute binary profile pattern.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
    
    padding : bool, default=False
        Pad sequences of unequal lengths with zeros at posterior end.

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 26*seq_length)
    
    Notes
    -----
    
    This algorithm is intended for proteins or peptides with equal lengths only.

    """

    # input handling
    X = check_input(X)

    # list of amino acids (IUPAC extended)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'

    # define maximum length 
    l = [len(seq) for seq in X]
    max_len = max(l)
    if padding == False and len(l) != l.count(max_len):
        raise ValueError('Sequences must be of equal length!')
        
    # compute binary profile pattern
    arr = np.zeros((len(X), len(amino_acids)*max_len))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        seq = seq[start-1:end] # positional information
        binary = [amino_acids.index(aa)+x*26 for x, aa in enumerate(seq)]
        arr[i,binary] = 1

    return arr