# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from .socn import socn
from ..utils.validation import check_input, check_alpha, check_natural

def qso(X, *, d=30, w=.1, remove_zero_cols=False, start=1, end=None): 
    """Quasi-sequence-order.

    This feature is derived from the distance matrix between the 20 amino acids.
    Here, we use two different distance matrices based on studies by Grantham 
    and Schneider-Wrede, respectively. For the exact method of calculating the
    quasi-sequence order, please refer to the documentation at
    https://protlearn.readthedocs.io/.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
        
    d : int, default=30
        Represents the lag. Must be smaller than sequence length.
        
    w : float, default=.1
        Weighting factor for the sequence-order effect.

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

    arr_sw :  ndarray of shape (n_samples, 20+d)
        Array containing QSO based on the Schneider-Wrede distance matrix.
    
    arr_g :  ndarray of shape (n_samples, 20+d)
        Array containing QSO based on the Grantham distance matrix.

    desc : list of length 20+d
        Order of QSO descriptors corresponding to columns in arr_sw and arr_g.
    
    References
    ----------

    Grantham, R. (1974). Amino acid difference formula to help explain 
    protein evolution. Science. 185 (4154): 862–864
    
    Schneider & Wrede (1994). The Rational Design of Amino Acid Sequences by 
    Artifical Neural Networks and Simulated Molecular Evolution: De Novo Design
    of an Idealized Leader Cleavge Site. Biophys Journal, 66, 335-344
    
    Chou,K.-C. (2000) Prediction of protein subcellar locations by incorporating
    quasi-sequence-order effect. Biochemical and Biophysical Research 
    Communications, 278, 477–483

    Examples
    --------

    >>> from protlearn.features import qso
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> sw, g, desc = qso(seqs, d=3, remove_zero_cols=True)
    >>> sw
    array([[0.66139001, 0.        , 0.        , 0.66139001, 0.66139001,
            0.        , 0.66139001, 0.66139001, 0.12285782, 0.15604737,
            0.0597048 ],
           [0.        , 1.33911348, 0.66955674, 0.66955674, 0.66955674,
            0.66955674, 0.66955674, 0.        , 0.07621852, 0.11022712,
            0.14399762]])
    >>> g
    array([[1.42887762e-04, 0.00000000e+00, 0.00000000e+00, 1.42887762e-04,
            1.42887762e-04, 0.00000000e+00, 1.42887762e-04, 1.42887762e-04,
            3.71008073e-01, 4.12445524e-01, 2.16403515e-01],
           [0.00000000e+00, 1.72010458e-04, 8.60052291e-05, 8.60052291e-05,
            8.60052291e-05, 8.60052291e-05, 8.60052291e-05, 0.00000000e+00,
            3.01095707e-01, 3.64610568e-01, 3.34207720e-01]])
    >>> desc
    ['A', 'E', 'G', 'K', 'L', 'P', 'R', 'Y', 'd1', 'd2', 'd3']

    """

    # input handling
    X = check_input(X)
    min_len = min([len(seq) for seq in X])
    if d >= min_len:
        raise ValueError('Lag parameter d must be smaller than sequence length!')
    
    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    desc = [aa for aa in amino_acids]
    
    for n in range(1, d+1):
        desc.append('d' + str(n))

    # calculate QSO
    arr_sw = np.zeros((len(X), len(amino_acids)+d))
    arr_g = np.zeros((len(X), len(amino_acids)+d))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
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
    if remove_zero_cols:
        cols_zeros = np.where(~arr_sw.any(axis=0))[0]
        arr_sw = np.delete(arr_sw, cols_zeros, axis=1)
        arr_g = np.delete(arr_g, cols_zeros, axis=1)
        desc = [i for j, i in enumerate(desc) if j not in cols_zeros]
        
    return arr_sw, arr_g, desc