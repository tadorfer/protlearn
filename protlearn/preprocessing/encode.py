# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from ..utils.validation import check_input

def encode(X, padding=False):
    """Label-encode amino acid sequences.

    The amino acids that serve as the building blocks for proteins and 
    peptides will be encoded into integers between 1-25. This is based on 
    the official IUPAC amino acid one-letter notation, which represents 22
    amino acids plus four additional symbols: B (asparagine or aspartic 
    acid), Z (glutamine or glutamic acid), J (leucine or isoleucine), and 
    X, which is used for unknown amino acids. Zeros are reserved for 
    optional padding. 
    
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
    amino_acids = IUPAC.ExtendedIUPACProtein().letters
    
    int_values = np.arange(1, len(amino_acids)+1)
    enc_list = []
    for seq in X:
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')
            
        seq = seq.upper()
        seq_split = [aa for aa in seq]
        seq_trans = [amino_acids.index(i)+1 for i in seq_split]
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