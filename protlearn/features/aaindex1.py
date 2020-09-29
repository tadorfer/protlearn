# Authors: Thomas Dorfer <thomas.a.dorfer@gmail.com>
#          Shoji Ihara <ihara@molcure.io>            

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def aaindex1(X, *, standardize='none', start=1, end=None):
    """AAIndex1-based physicochemical properties.

    AAindex1 ver.9.2 (release Feb, 2017) is a set of 20 numerical values
    representing various physicochemical and biological properties of amino
    acids. Currently, it contains 566 indices, of which 553 contain no NaNs.
    The indices will be collected for each amino acid in the sequence, 
    then averaged across the sequence. 

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.

    standardize : string, default='none'
        'none' : unstandardized index matrix will be returned
        'zscore' : index matrix is standardized to have
                   a mean of 0 and standard deviation of 1.
        'minmax' : index matrix is normalized to have a range of [0, 1].

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr : ndarray of shape (n_samples, 553-566) 
        Array containing the AAIndex1 physicochemical properties.

    desc : list of length 553-566
        Corresponds to the columns (AAIndices) in arr.

    Notes
    -----

    Columns (indices) containing NaNs will be removed. Thus, the resulting index
    matrix will have a column size between 553-566.

    References
    ----------

    Nakai, K., Kidera, A., and Kanehisa, M.; Cluster analysis of amino acid 
    indices for prediction of protein structure and function. Protein Eng. 2, 
    93-100 (1988). [PMID:3244698]

    Tomii, K. and Kanehisa, M.; Analysis of amino acid indices and mutation 
    matrices for sequence comparison and structure prediction of proteins. 
    Protein Eng. 9, 27-36 (1996). [PMID:9053899]

    Kawashima, S., Ogata, H., and Kanehisa, M.; AAindex: amino acid index 
    database. Nucleic Acids Res. 27, 368-369 (1999). [PMID:9847231]

    Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. 
    Nucleic Acids Res. 28, 374 (2000). [PMID:10592278]

    Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., 
    and Kanehisa, M.; AAindex: amino acid index database, progress report 2008. 
    Nucleic Acids Res. 36, D202-D205 (2008). [PMID:17998252]

    Examples
    --------

    >>> from protlearn.features import aaindex1
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> aaind, inds = aaindex1(seqs, standardize='zscore')
    >>> aaind.shape
    (2, 553)
    >>> len(inds)
    553

    """
    
    # input handling
    X = check_input(X)
    
    # Number of indices
    LEN = 566
    
    # list of amino acids (IUPAC extended)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    
    # load AAIndex1 data and get index names
    aaind1 = pd.read_csv(PATH+'aaindex1.csv')
    desc = aaind1['Description'].values
    
    # convert to dict for better performance
    aaind1 = {aa: aaind1[aa].to_numpy() for aa in amino_acids}

    # initialize empty array with shape (n_samples, 566)
    arr = np.zeros((len(X), LEN))

    # fill array with mean of indices per protein/peptide
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        tmp_arr = np.zeros((LEN, len(seq)) )

        # fill temporary array with indices for each amino acid of the sequence
        # and compute their mean across all rows
        for j, aa in enumerate(seq):
            tmp_arr[:,j] = aaind1[aa]

        # fill rows with mean vector of tmp_arr
        arr[i,:] = tmp_arr.mean(axis=1)

    # find columns with NaNs
    inds_nan = np.argwhere(np.isnan(arr))
    if len(inds_nan) != 0:
        cols_nan = np.unique(inds_nan[:,1])
        cols_nan = np.hstack(cols_nan)
    
    if standardize == 'none' or arr.shape[0] == 1:
        # remove columns with NaNs
        if len(inds_nan) != 0:
            arr = np.delete(arr, cols_nan, axis=1)
            desc = np.delete(desc, cols_nan)
        return arr, desc

    else:
        # remove columns with NaNs and all zeros
        cols_all = []
        cols_zeros = np.where(~arr.any(axis=0))[0]
            
        if len(inds_nan) != 0 and len(cols_zeros) != 0:
            cols_all = np.array((cols_nan, cols_zeros))
            cols_all = np.hstack(cols_all)

        elif len(inds_nan) != 0 and len(cols_zeros) == 0:
            cols_all = cols_nan

        elif len(inds_nan) == 0 and len(cols_zeros) != 0:
            cols_all = cols_zeros

        if len(cols_all) > 0:
            arr = np.delete(arr, cols_all, axis=1)
            desc = np.delete(desc, cols_all)

        # standardization
        if standardize == 'zscore':
            scaler = StandardScaler().fit(arr)
            arr = scaler.transform(arr)
            
            return arr, desc

        # normalization
        elif standardize == 'minmax':
            scaler = MinMaxScaler().fit(arr)
            arr = scaler.transform(arr)

            return arr, desc