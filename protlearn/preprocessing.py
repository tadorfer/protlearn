# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def check_input(X):
    """Check input format
    
    This function checks whether the input is in the correct format, i.e .either
    a string, a fasta object, or a list thereof. It also checks whether there are
    non-string types in the sequence and converts lower-case strings into upper-case
    strings.
    
    Parameters
    ----------
    
    X : string, fasta, or a list thereof
    
    Returns
    -------
    
    X_check : corrected data
        A ValueError is raised if incorrect data types are found.
    
    """ 
    
    ext = ['.fasta', '.faa', '.fa']
    err_msg = "Contains non-alphabetic characters"
    
    # list of strings
    if type(X) == list:
        print('list detected')
        X_check = []
        for seq in X:
            if str.isalpha(seq) == True:
                X_check.append(seq.upper())
            else:
                raise ValueError(err_msg)
                
    # strings: either sequence, of fasta file
    elif type(X) == str:
        _, extension = os.path.splitext(X)
        
        # fasta format
        if extension in ext:
            # single fasta sequence
            try:
                seq = SeqIO.read(X, 'fasta')
                seq = str(seq.seq)
                if str.isalpha(seq) == True:
                    X_check = seq.upper()
                else:
                    raise ValueError(err_msg)
                print('single fasta detected')
                
            # multiple fasta sequences
            except:
                data = [rec for rec in SeqIO.parse(X, 'fasta')]
                X_check = []
                for seq in data:
                    seq = str(seq.seq)
                    if str.isalpha(seq) == True:
                        X_check.append(seq.upper())
                    else:
                        raise ValueError(err_msg)
                print('multiple fasta detected')
        
        # string of sequence
        else:
            print('string detected')
            if str.isalpha(X) == True:
                X_check = X.upper()
            else:
                raise ValueError(err_msg)
            
    return X_check

def integer_encode(X, padding=False):
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

    Notes
    -----

    Amino acid sequence used for label-encoding were taken from the official
    IUPAC amino acid one-letter notation (Extended IUPAC Protein).

    """
    
    # list of amino acids for integer encoding (IUPAC extended)
    amino_acids = IUPAC.ExtendedIUPACProtein().letters
    
    ## input handling (single or multiple)
    # If input is a list of sequences (strings), nothing has to be specified.
    # If it is a string, it could either be a single sequence, or the path to
    # a fasta file. This block of code checks for that.
    ext = ['.fasta', '.faa', '.fa']
    if type(X) == str:
        _, extension = os.path.splitext(X)
        
        # fasta format
        if extension in ext:
            # single fasta sequence
            try:
                X = [str(SeqIO.read(X, 'fasta').seq)]
                
            # multiple fasta sequences
            except:
                X = [str(rec.seq) for rec in SeqIO.parse(X, 'fasta')]
                
        else:
            X = [X] 
    
    int_values = np.arange(1, len(amino_acids)+1)
    enc_list = []
    for seq in X:
        # format must be string
        if type(seq) != str:
            seq = str(seq.seq)
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data type must be string!')
            
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
        return enc_arr[0]
    else:
        return enc_arr