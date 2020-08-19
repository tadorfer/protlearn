# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from Bio.Alphabet import IUPAC
from utils.validation import check_input

def binary(X, start=1, end=None):
    """Compute binary profile pattern.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

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
    amino_acids = IUPAC.ExtendedIUPACProtein.letters

    # compute binary profile pattern
    arr = np.zeros((len(X), len(amino_acids)*len(X[0])))
    for i, seq in enumerate(X):
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')

        seq = seq[start-1:end] # positional information
        binary = [amino_acids.index(aa)+x*26 for x, aa in enumerate(seq)]
        arr[i,binary] = 1

    return arr