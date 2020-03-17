# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pkg_resources


PATH = pkg_resources.resource_filename('protlearn', 'docs/')

def length(X, method='int'):
    """Compute the length of proteins or peptides.
    
    The number of amino acids that a protein or peptide is comprised of will be
    counted and returned as either an integer array or a one-hot-encoded array.

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    method : string, default='int'

        'int' : compute array of protein/peptide lengths
        'ohe' : compute one-hot encoded array of protein/peptide lengths

    Returns
    -------

    lengths : ndarray of shape (n_samples, ) if method = 'int'.
              ndarray of shape (n_samples, n_unique_lengths) if method = 'ohe'

    """

    # list of protein/peptide lengths
    all_len = [len(X['Sequence'][i]) for i in range(X.shape[0])]
    
    # one-hot-encoded lengths of proteins/peptides
    len_span = max(all_len) - min(all_len)
    len_ohe = np.zeros((X.shape[0], len_span+1))
    len_unique = np.arange(min(all_len), max(all_len)+1)

    # fill array with ones based on sequence lengths (columns)
    for i in range(len(all_len)):    
        seq_len = all_len[i]
        j = np.where(len_unique == seq_len)[0][0]
        len_ohe[i,j] = 1

    # delete columns with all zeros
    zero_cols = np.where(~len_ohe.any(axis=0))[0]
    len_ohe = np.delete(len_ohe, zero_cols, axis=1)

    return {
            'int': lambda: np.asarray(all_len), 
            'ohe': lambda: len_ohe,
    }.get(method, lambda: None)()


def composition(X, method='absolute', round_fraction=3):
    """Compute the amino acid composition of proteins or peptides.

    The frequency of each of the 20 amino acids in a protein or peptide are 
    counted and returned using the following column order:

    'A','C','D','E','F','G','H','I','K','L',
    'M','N','P','Q','R','S','T','V','W','Y'

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    method : string, default='absolute'

        'absolute' : compute absolute amino acid composition
        'relative' : compute relative amino acid composition

    round_fraction : int, default=3
        This applies only if method='relative'. For shorter peptides with only 
        a few amino acids, the rounding of amino acid fractions will not have a 
        big impact. However, with longer proteins, decimal places become 
        increasingly more significant and can thus be increased.

    Returns
    -------

    comp :  Pandas DataFrame of shape (n_samples, n_unique_amino_acids)

    Notes
    -----

    The returned dataframe 'comp' can easily be converted into a numpy array 
    with the command 'np.asarray(comp)', if desired.

    """

    # list of amino acids repesenting array columns
    amino_acids = ['A','C','D','E','F','G','H','I','K','L',
                   'M','N','P','Q','R','S','T','V','W','Y']

    # initialize empty array with shape (n_samples, 20)
    comp_arr = np.zeros((X.shape[0], len(amino_acids)))

    # fill array with frequency of amino acids per protein/peptide
    for i in range(X.shape[0]):
        seq = [aa for aa in X['Sequence'][i]]
        values, counts = np.unique(seq, return_counts=True)
        indices = [amino_acids.index(j) for j in values]
        comp_arr[i,indices] = counts

    # delete zero columns
    cols_zeros = np.where(~comp_arr.any(axis=0))[0]
    comp_arr = np.delete(comp_arr, cols_zeros, axis=1)
    amino_acids = [i for j, i in enumerate(amino_acids) if j not in cols_zeros]

    if method == 'absolute':
        return pd.DataFrame(comp_arr, columns=amino_acids)

    elif method == 'relative':
        all_len = [len(X['Sequence'][i]) for i in range(X.shape[0])]
        comp_rel = [np.round((comp_arr[i]/all_len[i]), round_fraction)\
                    for i in range(comp_arr.shape[0])]
        
        return pd.DataFrame(np.asarray(comp_rel), columns=amino_acids)
        

def aaindex1(X, standardize='none'):
    """Compute amino acid indices from AAIndex1.

    AAindex1 ver.9.2 (release Feb, 2017) is a set of 20 numerical values
    representing different physicochemical and biological properties of amino
    acids. Currently, it contains 566 such indices, of which 553 contain no 
    NaNs. The indices will be collected for each amino acid in the sequence, 
    then averaged across the sequence. 

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    standardize : string, default='none'

        'none' : unstandardized index matrix will be returned
        'zscore' : index matrix is standardized across columns (indices) to have
                   a mean of 0 and standard deviation of 1 (unit variance).
        'minmax' : index matrix is scaled (normalized) across columns (indices)
                   to have a range of [0, 1].

    Returns
    -------

    arr_index1 : ndarray of shape (n_samples, 553) 

    Notes
    -----

    Columns (indices) containing NaNs will be removed. Thus, the resulting index
    matrix will have a column size of 553, rather than 566.

    The returned dataframe 'arr_index1' can easily be converted into a numpy 
    array with the command 'np.asarray(arr_index1)', if desired.

    """
    
    # load AAIndex1 data
    aaind1 = pd.read_csv(PATH+'aaindex1.csv')

    # get descriptions of all 566 indices
    desc = aaind1['Description'].values

    # initialize empty array with shape (n_samples, n_aa1_indices)
    aaind_arr = np.zeros((len(X), aaind1.shape[0]))

    # fill array with mean of indices per protein/peptide
    for i, seq in enumerate(range(len(X))):
        sequence = X['Sequence'][seq]
        tmp_arr = np.zeros((aaind1.shape[0], len(sequence)) )

        # fill temporary array with indices for each amino acid of the sequence
        # and compute their mean across all rows
        for j, aa in enumerate(sequence):
            tmp_arr[:,j] = np.asarray(aaind1[aa])

        # fill rows with mean vector of tmp_arr
        aaind_arr[i,:] = tmp_arr.mean(axis=1)

    if standardize == 'none':
        # removing columns with NaNs
        inds_nan = np.argwhere(np.isnan(aaind_arr))
        if len(inds_nan) != 0:
            cols_nan = np.unique(inds_nan[:,1])
            cols_nan = np.hstack(cols_nan)
            aaind_arr = np.delete(aaind_arr, cols_nan, axis=1)
            desc = np.delete(desc, cols_nan)
        return pd.DataFrame(aaind_arr, columns=desc)

    else:
        # finding and removing columns with NaNs and all zeros
        cols_all = []
        inds_nan = np.argwhere(np.isnan(aaind_arr))
        if len(inds_nan) != 0:
            cols_nan = np.unique(inds_nan[:,1])
            cols_nan = np.hstack(cols_nan)

        cols_zeros = np.where(~aaind_arr.any(axis=0))[0]
            
        if len(inds_nan) != 0 and len(cols_zeros) != 0:
            cols_all = np.array((cols_nan, cols_zeros))
            cols_all = np.hstack(cols_all)

        elif len(inds_nan) != 0 and len(cols_zeros) == 0:
            cols_all = cols_nan

        elif len(inds_nan) == 0 and len(cols_zeros) != 0:
            cols_all = cols_zeros

        if len(cols_all) > 0:
            aaind_arr = np.delete(aaind_arr, cols_all, axis=1)
            desc = np.delete(desc, cols_all)

        # standardization
        if standardize == 'zscore':
            scaler = StandardScaler().fit(aaind_arr)
            aaind_arr = scaler.transform(aaind_arr)
            
            return pd.DataFrame(aaind_arr, columns=desc)

        # normalization
        elif standardize == 'minmax':
            scaler = MinMaxScaler().fit(aaind_arr)
            aaind_arr = scaler.transform(aaind_arr)

            return pd.DataFrame(aaind_arr, columns=desc)
        
        
def aaindex2(X, standardize='none'):
    """Compute amino acid indices from AAIndex2.

    AAindex2 ver.9.2 (release Feb, 2017) is a set of lower triangular (67), 
    square (21), and rectangular (6) matrices and includes a collection of 
    published amino acid substitution matrices. Currently, it contains 
    94 such matrices. Some of the square and rectangular matrices include
    gaps and various metrics for cysteines (disulfide-bonded and free), which
    are removed for all computations (only the standard 20 amino acids are
    used here). Also, two indices contain NaNs. Therefore, if any of the 
    sequences contains an amino acid pair whose index is NaN, the entire
    column will be removed.
    
    The substitution score will be collected for each amino acid pair in the 
    sequence, then averaged across the sequence.

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    standardize : string, default='none'

        'none' : unstandardized index matrix will be returned
        'zscore' : index matrix is standardized across columns (indices) to have
                   a mean of 0 and standard deviation of 1 (unit variance).
        'minmax' : index matrix is scaled (normalized) across columns (indices)
                   to have a range of [0, 1].

    Returns
    -------

    arr_index2 : ndarray of shape (n_samples, 92-94) if standardize='none'
        Column size could vary when standardize != 'none'.

    Notes
    -----

    Columns (indices) containing NaNs will be removed. Thus, the resulting index
    matrix will have a column size between 92-94.

    The returned dataframe 'arr_index2' can easily be converted into a numpy 
    array with the command 'np.asarray(arr_index1)', if desired.

    """
    
    # load AAIndex2 data
    raw_lt = pd.read_csv(PATH+'aaindex2_lowtri.csv')
    raw_sq = pd.read_csv(PATH+'aaindex2_square.csv')

    def compute_aaind2(index, shape):
        "Compute AAIndex2 for all shapes"
        
        # get unique index descriptions
        desc = [index['Description'][i] for i in\
                np.arange(0, index.shape[0], 20)]
        
        # drop columns to facilitate easier indexing
        index = index.drop(['Description', 'Amino Acids'], 
                            axis=1).reset_index(drop=True)
        
        # get lengths and indices (e.g. 0, 20, 40,...)
        LENGTH = len(desc)
        inds = np.arange(0, index.shape[0], 20)
        
        # compute dict for indexing ({'A': 0, 'R': 1, ...})
        aa_dict = {index.columns[i]: i for i in range(index.shape[1])}
        
        # initialitze empty array
        arr = np.zeros((len(X), LENGTH))
        
        for a, seq in enumerate(range(len(X))):
            sequence = X['Sequence'][seq]
            conpot = np.zeros((len(sequence)-1, LENGTH))

            for i in range(len(sequence)-1):
                aa1, aa2 = sequence[i], sequence[i+1]
                aa1_ind, aa2_ind = aa_dict[aa1], aa_dict[aa2]

                if shape == 'square':
                    aa1_ind = aa1_ind + inds
                    conpot[i,:] = index.iloc[aa1_ind, aa2_ind]

                elif shape == 'lowtri':
                    if aa1_ind > aa2_ind:
                        aa1_ind += inds
                        conpot[i,:] = index.iloc[aa1_ind, aa2_ind]
                    else:
                        aa2_ind += inds
                        conpot[i,:] = index.iloc[aa2_ind, aa1_ind]

            arr[a,:] = np.asarray(conpot.mean(axis=0))
            
        if standardize == 'none':
            # removing columns with NaNs
            inds_nan = np.argwhere(np.isnan(arr))
            if len(inds_nan) != 0:
                cols_nan = np.unique(inds_nan[:,1])
                cols_nan = np.hstack(cols_nan)
                arr = np.delete(arr, cols_nan, axis=1)
                desc = np.delete(desc, cols_nan)
            return pd.DataFrame(arr, columns=desc)
        
        else:
            # finding and removing columns with NaNs and all zeros
            cols_all = []
            inds_nan = np.argwhere(np.isnan(arr))
            if len(inds_nan) != 0:
                cols_nan = np.unique(inds_nan[:,1])
                cols_nan = np.hstack(cols_nan)

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
                
                return pd.DataFrame(arr, columns=desc)
            
            elif standardize == 'minmax':
                scaler = MinMaxScaler().fit(arr)
                arr = scaler.transform(arr)
                
                return pd.DataFrame(arr, columns=desc)
            
        
    # get aaindices2
    lt = compute_aaind2(raw_lt, 'lowtri')
    sq = compute_aaind2(raw_sq, 'square')
    
    return pd.concat([lt, sq], axis=1)


def aaindex3(X, standardize='none'):
    """Compute amino acid indices from AAIndex3.

    AAindex3 ver.9.2 (release Feb, 2017) is a set of lower triangular (44)
    and square (3) matrices and includes a collection of published protein 
    pairwise contact potentials. Currently, it contains 47 such matrices, of 
    which four contain NaNs. Therefore, if any of the sequences contains an 
    amino acid pair whose index is NaN, the entire column will be removed.

    The pairwise contact potentials will be collected for each amino acid pair
    in the sequence, then averaged across the sequence.

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    standardize : string, default='none'

        'none' : unstandardized index matrix will be returned
        'zscore' : index matrix is standardized across columns (indices) to have
                   a mean of 0 and standard deviation of 1 (unit variance).
        'minmax' : index matrix is scaled (normalized) across columns (indices)
                   to have a range of [0, 1].

    Returns
    -------

    arr_index3 : ndarray of shape (n_samples, 43-47) if standardize='none'
        Column size could vary when standardize != 'none'.

    Notes
    -----

    Columns (indices) containing NaNs will be removed. Thus, the resulting index
    will have a column size between 43-47. 

    The returned dataframe 'arr_index3' can easily be converted into a numpy 
    array with the command 'np.asarray(arr_index3)', if desired.

    """
    
    # load AAIndex3 data
    raw_lt = pd.read_csv(PATH+'aaindex3_lowtri.csv')
    raw_sq = pd.read_csv(PATH+'aaindex3_square.csv')

    def compute_aaind3(index, shape):
        "Compute AAIndex3 for all shapes"
        
        # get unique index descriptions
        desc = [index['Description'][i] for i in\
                      np.arange(0, index.shape[0], 20)]
        
        # drop columns to facilitate easier indexing
        index = index.drop(['Description', 'Amino Acids'], 
                            axis=1).reset_index(drop=True)
        
        # get lengths and indices (e.g. 0, 20, 40,...)
        LENGTH = len(desc)
        inds = np.arange(0, index.shape[0], 20)
        
        # compute dict for indexing ({'A': 0, 'R': 1, ...})
        aa_dict = {index.columns[i]: i for i in range(index.shape[1])}
        
        # initialitze empty array
        arr = np.zeros((len(X), LENGTH))
        
        for a, seq in enumerate(range(len(X))):
            sequence = X['Sequence'][seq]
            conpot = np.zeros((len(sequence)-1, LENGTH))

            for i in range(len(sequence)-1):
                aa1, aa2 = sequence[i], sequence[i+1]
                aa1_ind, aa2_ind = aa_dict[aa1], aa_dict[aa2]

                if shape == 'square':
                    aa1_ind = aa1_ind + inds
                    conpot[i,:] = index.iloc[aa1_ind, aa2_ind]

                elif shape == 'lowtri':
                    if aa1_ind > aa2_ind:
                        aa1_ind += inds
                        conpot[i,:] = index.iloc[aa1_ind, aa2_ind]
                    else:
                        aa2_ind += inds
                        conpot[i,:] = index.iloc[aa2_ind, aa1_ind]

            arr[a,:] = np.asarray(conpot.mean(axis=0))
            
        if standardize == 'none':
            # removing columns with NaNs
            inds_nan = np.argwhere(np.isnan(arr))
            if len(inds_nan) != 0:
                cols_nan = np.unique(inds_nan[:,1])
                cols_nan = np.hstack(cols_nan)
                arr = np.delete(arr, cols_nan, axis=1)
                desc = np.delete(desc, cols_nan)
            return pd.DataFrame(arr, columns=desc)
        
        else:
            # finding and removing columns with NaNs and all zeros
            cols_all = []
            inds_nan = np.argwhere(np.isnan(arr))
            if len(inds_nan) != 0:
                cols_nan = np.unique(inds_nan[:,1])
                cols_nan = np.hstack(cols_nan)

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
                
                return pd.DataFrame(arr, columns=desc)
            
            elif standardize == 'minmax':
                scaler = MinMaxScaler().fit(arr)
                arr = scaler.transform(arr)
                
                return pd.DataFrame(arr, columns=desc)
            
        
    # get aaindices3
    lt = compute_aaind3(raw_lt, 'lowtri')
    sq = compute_aaind3(raw_sq, 'square')
    
    return pd.concat([lt, sq], axis=1)


def ngram_composition(X, ngram=2):
    """Compute n-gram peptide composition.
    
    This function computes the di-, tri-, or quadpeptide composition of
    amino acid sequences. Therefore, the function argument 'ngram' can only
    take on the values 2, 3, and 4 - otherwise, it will raise a ValueError.
    
    Parameters
    ----------
    
    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.
       
    ngram : int, default=2
        Integer denoting the desired n-gram composition.
        
        2 : dipeptide composition
        3 : tripepitde composition
        4 : quadpeptide composition
        
    Returns
    -------
    
    arr_ngram : ndarray of shape (n_samples, n_unique_20^ngram)
        Depending on ngram, returned array will be of size:
        - (n_samples, 400) for dipeptide composition
        - (n_samples, 8000) for tripeptide composition
        - (n_samples, 160000) for quadpeptide composition

    Notes
    -----

    Columns containing all zeros will be removed, therefore the column size of
    the returned array can vary.
       
    """
    
    # make sure ngram is between 2-4
    valid = [2, 3, 4]
    if ngram not in valid:
        raise ValueError("ngram_comp: ngram must be one of %r." % valid)
    
    # get ngram combinations
    aa_combo = []
    amino_acids = ['A','C','D','E','F','G','H','I','K','L',
                   'M','N','P','Q','R','S','T','V','W','Y']

    def combo(seq, prefix, n, k):
        if k == 0: 
            aa_combo.append(prefix)
            return aa_combo  

        for i in range(n): 
            newPrefix = prefix + seq[i] 
            combo(seq, newPrefix, n, k-1) 

    combo(amino_acids, "", 20, ngram)
    
    # create dataframe with all zeros
    arr_ngram = np.zeros((len(X), len(aa_combo)))
    arr_ngram = pd.DataFrame(arr_ngram, columns=aa_combo)
    
    # define n-gram function
    def n_gram(seq):
        "Computing n-gram features for sequence pairs"
        n_gram_list = []
        for i in range(len(seq)):
            aa_pair = seq[i:i+ngram]
            if len(aa_pair) == ngram:
                n_gram_list.append(aa_pair)
        return Counter(n_gram_list)
    
    # compute n-gram composition
    for i in range(len(X)):
        pep_comp = n_gram(X['Sequence'][i])
        for j in range(len(pep_comp)):
            keys = list(pep_comp)
            arr_ngram[keys[j]][i] = pep_comp[keys[j]]

    # delete zero columns
    arr_ngram = arr_ngram.loc[:, (arr_ngram!=0).any(axis=0)]
            
    return arr_ngram