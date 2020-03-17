# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd


def txt_to_df(X, label=None):
    """Convert .txt file containing peptide sequences into Pandas DataFrame. 

    The .txt file must contain one peptide sequence per row and no header.
    For classification tasks, the second function argument 'label' takes in an
    integer class label and generates a column with this class label.

    Parameters
    ----------

    X : .txt file 
        Contains peptide sequences (one per row) and no header.
    
    label : int, default=None
        Integer label denoting class (optional). 

    Returns
    -------

    df: DataFrame of shape (n_samples, 2) with columns 'Sequence' and 'Label'

    """

    df = pd.read_csv(X, names=['Sequence'], header=None)

    if label:
        df['Label'] = np.ones(df.shape[0])*label
        
    return df


def integer_encode(X, padding=False):
    """Label-encode amino acid sequences.

    The 20 amino acids that serve as the building blocks for proteins and 
    peptides will be encoded into integers between 1-20. Zeros are reserved
    for optional padding. Integer-encoding can be useful for the classification
    of proteins or peptides using sequence-based prediction models such as 
    LSTMs or GRUs. 

    Parameters
    ----------

    X : Pandas DataFrame
        Columns must contain 'Sequence' and 'Label'.

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

    Amino acid sequence used for label-encoding:

    'A','C','D','E','F','G','H','I','K','L',
    'M','N','P','Q','R','S','T','V','W','Y'

    """

    # list of amino acids for integer encoding
    amino_acids = ['A','C','D','E','F','G','H','I','K','L',
                   'M','N','P','Q','R','S','T','V','W','Y'] 
    int_values = np.arange(1, len(amino_acids)+1)

    enc_list = []
    for seq in X['Sequence']:
        seq_split = [aa for aa in seq]
        seq_trans = [amino_acids.index(i)+1 for i in seq_split]
        enc_list.append(np.asarray(seq_trans))
    enc_arr = np.asarray(enc_list)

    if padding == True:
        # get maximum sequence length
        all_len = [len(i) for i in enc_arr]
        max_len = max(all_len)

        # padding up to max_len with zeros
        enc_arr = [np.pad(enc_arr[i], (0, max_len-len(enc_arr[i])),\
                  'constant') for i in range(enc_arr.shape[0])]
        enc_arr = np.asarray(enc_arr)

    return enc_arr