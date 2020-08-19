# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import re
import numpy as np
from itertools import product
from Bio.Alphabet import IUPAC
from utils.validation import check_input

def cksaap(X, k=1, start=1, end=None):
    """Compute composition of k-spaced amino acid pairs.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
    
    lambda_ : int, default=1
        Counted rank (tier) of the correlation along an amino acid sequence.
        
    k : int, default=1
        Space between two amino acid pairs.

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 676)
    
    patterns : amino acid pairs with k gaps

    """
    
    # input handling
    X = check_input(X)
    
    # initialize empty array 
    arr = np.empty((len(X), 676), dtype=int)
    
    # list of amino acids (IUPAC extended)
    amino_acids = IUPAC.ExtendedIUPACProtein.letters
    doublets = sorted([c[0]+c[1] for c in product(amino_acids, repeat=2)])
    patterns = [doublets[i][0]+'.'*k+doublets[i][1] for i in range(len(doublets))]

    # compute CKSAAP
    for i, seq in enumerate(X):
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')
            
        seq = seq[start-1:end] # positional information
        for j, pattern in enumerate(patterns):
            cnt_pattern = len(re.findall(r'(?=('+pattern+'))', seq))
            arr[i, j] = cnt_pattern
            
    # delete zero columns
    cols_zeros = np.where(~arr.any(axis=0))[0]
    arr = np.delete(arr, cols_zeros, axis=1)
    patterns = [i for j, i in enumerate(patterns) if j not in cols_zeros]
    
    return arr, patterns