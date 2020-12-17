# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import re
import numpy as np
from itertools import product
from ..utils.validation import check_input, check_alpha, check_natural

def cksaap(X, *, k=1, remove_zero_cols=False, start=1, end=None):
    """Composition of k-spaced amino acid pairs.

    This function returns the k-spaced amino acid pair composition of each 
    sequence in the dataset. Since there are 20 natural amino acids, there are 
    400 possible amino acid pairs. The parameter 'k' represents the gap between
    the amino acid pair. An example for k=1 would be AxY, where 'x' can be any 
    amino acid. Similary, an example for k=2 would be AxxY. If k=0, the function
    returns the dipeptide composition of each sequence.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
    
    lambda_ : int, default=1
        Counted rank (tier) of the correlation along an amino acid sequence.
        
    k : int, default=1
        Space between two amino acid pairs.

    remove_zero_cols : bool, default=False
        If true, columns containing only zeros will be deleted. 

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 400)
        Array containing k-spaced amino acid pair composition.
    
    patterns : list of length 400
        Amino acid pairs with k gaps corresponding to columns in arr.

    References
    ----------

    Chen, K., Kurgan, L.A. & Ruan, J. Prediction of flexible/rigid regions from
    protein sequences using k-spaced amino acid pairs. BMC Struct Biol 7, 25 
    (2007). https://doi.org/10.1186/1472-6807-7-25

    Examples
    --------

    >>> from protlearn.features import cksaap
    >>> seqs = ['ARKLY', 'EERKPGL', 'AAAAAALY']
    >>> ck, pairs = cksaap(seqs, remove_zero_cols=True)
    >>> ck
    array([[0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0],
           [0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1],
           [4, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]])
    >>> pairs
    ['A.A', 'A.K', 'A.L', 'A.Y', 'E.K', 'E.R', 'K.G', 'K.Y', 'P.L', 'R.L', 'R.P']
    >>> ck2, pairs2 = cksaap(seqs, k=2, remove_zero_cols=True)
    >>> ck2
    array([[0, 1, 0, 0, 0, 0, 0, 1],
           [0, 0, 0, 1, 1, 1, 1, 0],
           [3, 1, 1, 0, 0, 0, 0, 0]])
    >>> pairs2
    ['A..A', 'A..L', 'A..Y', 'E..K', 'E..P', 'K..L', 'R..G', 'R..Y']

    """
    
    # input handling
    X = check_input(X)
    
    # initialize empty array 
    arr = np.empty((len(X), 400), dtype=int)
    
    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    doublets = sorted([c[0]+c[1] for c in product(amino_acids, repeat=2)])
    patterns = [doublets[i][0]+'.'*k+doublets[i][1] for i in range(len(doublets))]

    # compute CKSAAP
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        for j, pattern in enumerate(patterns):
            cnt_pattern = len(re.findall(r'(?=('+pattern+'))', seq))
            arr[i, j] = cnt_pattern
            
    # delete zero columns
    if remove_zero_cols:
        cols_zeros = np.where(~arr.any(axis=0))[0]
        arr = np.delete(arr, cols_zeros, axis=1)
        patterns = [i for j, i in enumerate(patterns) if j not in cols_zeros]
    
    return arr, patterns