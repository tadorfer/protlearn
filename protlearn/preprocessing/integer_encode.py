# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha

def integer_encode(X, notation='standard', padding=False):
    """Encode amino acids as integers.

    This function converts amino acids into their corresponding integers 
    based on the specified notation, starting at 1. Zeros are reserverd for 
    optional padding. This is particularly useful for preparing a sequence-based
    model such as a long short-term memory (LSTM) or a gated recurrent unit (GRU). 

    Parameters
    ----------

    X : string, fasta, or a list thereof
        Dataset of amino acid sequences.

    notation : string, default='standard'
        'standard' : 20 natural amino acids
        'extended' : 20 natural + six additional amino acids:
            B = aspartic acid or asparagine
            X = unknown or 'other' amino acid
            Z = glutamic acid or glutamine
            J = leucine or isoleucine
            U = selenocysteine
            O = pyrrolysine

    padding : bool, default=False
        False : sequences are returned in their original lengths
        True : sequences will be padded with zeros at the end up until the
               length of the longest sequence in the dataset

    Returns
    -------

    enc : ndarray of shape (n_samples,) if padding=False
          ndarray of shape (n_samples, max_len) if padding=True
        Contains the label-encoded peptide sequences.

    amino_acids : amino acid order of enc array
        This serves as a lookup for the encoded sequences.

    Examples
    --------

    >>> from protlearn.preprocessing import integer_encode
    >>> seqs = ['ARKLY', 'EERNPAA', 'QEPGPGLLLK']
    >>> enc, aa = integer_encode(seqs, padding=True)
    >>> enc
    array([[ 1, 15,  9, 10, 20,  0,  0,  0,  0,  0],
           [ 4,  4, 15, 12, 13,  1,  1,  0,  0,  0],
           [14,  4, 13,  6, 13,  6, 10, 10, 10,  9]])
    >>> aa
    'ACDEFGHIKLMNPQRSTVWY'

    Notes
    -----

    Amino acid sequence used for label-encoding were taken from the official
    IUPAC amino acid one-letter notation (Extended IUPAC Protein).

    """
    
    # input handling 
    X = check_input(X)

    # list of amino acids for integer encoding (IUPAC extended)
    if notation == 'standard':
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    elif notation == 'extended':
        amino_acids = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
    else:
        raise ValueError("Valid options are 'standard' and 'extended'.")
    
    aa_dict = {amino_acids[i]: i+1 for i in range(len(amino_acids))}
    enc_list = []
    for seq in X:
        check_alpha(seq) # check if alphabetical      
        seq = seq.upper()
        seq_trans = [aa_dict[aa] for aa in seq]
        enc_list.append(np.asarray(seq_trans))
    enc_arr = np.asarray(enc_list)

    if padding == True and enc_arr.shape[0] > 1:
        # get maximum sequence length
        all_len = [len(i) for i in enc_arr]
        max_len = max(all_len)

        # padding up to max_len with zeros
        enc_arr = [np.pad(enc_arr[i], (0, max_len-len(enc_arr[i])),\
                  'constant') for i in range(enc_arr.shape[0])]
        enc_arr = np.asarray(enc_arr)

    if enc_arr.shape[0] == 1:
        return enc_arr[0], amino_acids
    else:
        return enc_arr, amino_acids