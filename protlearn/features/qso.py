# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from Bio.Alphabet import IUPAC
from .socn import socn
from ..utils.validation import check_input, check_alpha

def qso(X, d=30, w=.1, start=1, end=None): 
    """Compute Quasi-sequence-order.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        
    d : int, default=30
        Represents the lag. Must be smaller than sequence length.
        
    w : float, default=.05
        Weight factor.
    
    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr_sw :  ndarray of shape (n_samples, 20+d)
    
    arr_g :  ndarray of shape (n_samples, 20+d)

    """

    # input handling
    X = check_input(X)
    min_len = min([len(seq) for seq in X])
    if d >= min_len:
        raise ValueError('Lag parameter d must be smaller than sequence length!')
    
    # list of amino acids (IUPAC standard)
    amino_acids = IUPAC.IUPACProtein.letters
    desc = [aa for aa in amino_acids]
    
    for n in range(1, d+1):
        desc.append('d' + str(n))

    # calculate QSO
    arr_sw = np.zeros((len(X), len(amino_acids)+d))
    arr_g = np.zeros((len(X), len(amino_acids)+d))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        seq = seq[start-1:end] # positional information 
        socn_sw = socn(seq, d=d)[0]
        socn_g = socn(seq, d=d)[1]
        
        cnt = Counter(seq)
        
        qso_sw = [cnt[aa] / (1+w * sum(socn_sw)) for aa in amino_acids]
        qso_sw = qso_sw + [(w*j) / (1+w * sum(socn_sw)) for j in socn_sw]
        arr_sw[i,:] = qso_sw
        
        qso_g = [cnt[aa] / (1+w * sum(socn_g)) for aa in amino_acids]
        qso_g = qso_g + [(w*j) / (1+w * sum(socn_g)) for j in socn_g]
        arr_g[i,:] = qso_g

    # delete zero columns
    cols_zeros = np.where(~arr_sw.any(axis=0))[0]
    arr_sw = np.delete(arr_sw, cols_zeros, axis=1)
    arr_g = np.delete(arr_g, cols_zeros, axis=1)
    desc = [i for j, i in enumerate(desc) if j not in cols_zeros]
        
    return arr_sw, arr_g, desc