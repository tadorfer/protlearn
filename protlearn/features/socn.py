# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def socn(X, *, d=30, start=1, end=None): 
    """Sequence-order-coupling number.

    This feature is derived from the distance matrix between the 20 amino acids.
    Here, we use two different distance matrices based on studies by Grantham 
    and Schneider-Wrede, respectively. For the exact method of calculating the
    sequence-order-coupling number, please refer to the documentation at
    https://protlearn.readthedocs.io/.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
        
    d : int, default=30
        Represents the lag. Must be smaller than sequence length.
    
    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr_sw :  ndarray of shape (n_samples, d)
        Array containing SOCN based on the Schneider-Wrede distance matrix.
    
    arr_g :  ndarray of shape (n_samples, d)
        Array containing SOCN based on the Grantham distance matrix.

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

    >>> from protlearn.features import socn
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> sw, g = socn(seqs, d=3)
    >>> sw
    array([1.85757 , 2.359385, 0.902717],
          [1.138343, 1.64627 , 2.150641]])
    >>> g
    array([[25965., 28865., 15145.],
           [35009., 42394., 38859.]])

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
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_dict = {amino_acids[i]: i for i in range(20)}

    # calculate SOCN
    arr_sw = np.zeros((len(X), d))
    arr_g = np.zeros((len(X), d))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical
        check_natural(seq) # check for unnatural amino acids  
        seq = seq[start-1:end] # positional information 
        tmp_sw = np.zeros((d,))
        tmp_g = np.zeros((d,))
        for n, lag in enumerate(range(1, d+1)):
            tmp_sw[n] = sum([(sw[aa_dict[seq[j]], aa_dict[seq[j+lag]]])**2 \
                            for j in range(len(seq)-lag)])
            tmp_g[n] = sum([(g[aa_dict[seq[j]], aa_dict[seq[j+lag]]])**2 \
                            for j in range(len(seq)-lag)])
        
        arr_sw[i,:] = tmp_sw
        arr_g[i,:] = tmp_g
        
    return arr_sw, arr_g