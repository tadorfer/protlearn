# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha

def integer_encode(X, padding=False, notation='standard'):
    """Label-encode amino acid sequences.

    The amino acids  will be encoded into integers given the chosen notation.
    Zeros are reserved for optional padding. 
    
    Integer-encoding can be useful for the classification of proteins or 
    peptides using sequence-based prediction models such as LSTMs or GRUs. 

    Parameters
    ----------

    X : string, fasta, or a list thereof

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