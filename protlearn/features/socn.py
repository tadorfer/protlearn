# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from Bio.Alphabet import IUPAC
from ..utils.validation import check_input
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def socn(X, d=30, start=1, end=None): 
    """Compute Sequence-order-coupling number.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        
    d : int, default=30
        Represents the lag. Must be smaller than sequence length.
    
    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr_sw :  ndarray of shape (n_samples, d)
    
    arr_g :  ndarray of shape (n_samples, d)

    """

    # input handling
    X = check_input(X)
    min_len = min([len(seq) for seq in X])
    if d >= min_len:
        raise ValueError('Lag parameter d must be smaller than sequence length!')

    # load data
    df_sw = pd.read_csv(PATH+'schneider-wrede.csv').set_index('AminoAcid')
    df_g = pd.read_csv(PATH+'grantham.csv').set_index('AminoAcid')
    sw = np.asarray(df_sw)
    g = np.asarray(df_g)
    
    # list of amino acids (IUPAC standard)
    amino_acids = IUPAC.IUPACProtein.letters
    aadict = {amino_acids[i]: i for i in range(20)}

    # calculate SOCN
    arr_sw = np.zeros((len(X), d))
    arr_g = np.zeros((len(X), d))
    for i, seq in enumerate(X):
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')

        seq = seq[start-1:end] # positional information 
        tmp_sw = np.zeros((d,))
        tmp_g = np.zeros((d,))
        for n, lag in enumerate(range(1, d+1)):
            tmp_sw[n] = sum([(sw[aadict[seq[j]], aadict[seq[j+lag]]])**2 for j in range(len(seq)-lag)])
            tmp_g[n] = sum([(g[aadict[seq[j]], aadict[seq[j+lag]]])**2 for j in range(len(seq)-lag)])
        
        arr_sw = tmp_sw
        arr_g = tmp_g
        
    return arr_sw, arr_g