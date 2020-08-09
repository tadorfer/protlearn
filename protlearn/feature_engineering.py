# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
import pandas as pd
from collections import Counter
from Bio.Alphabet import IUPAC
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pkg_resources


PATH = pkg_resources.resource_filename('protlearn', 'docs/')

def length(X, method='int'):
    """Compute the length of proteins or peptides.
    
    The number of amino acids that a protein or peptide is comprised of will be
    counted and returned as either an integer (single sequence), array of integers,
    or a one-hot-encoded array.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    method : string, default='int'

        'int' : compute array of protein/peptide lengths
        'ohe' : compute one-hot encoded array of protein/peptide lengths

    Returns
    -------

    lengths : ndarray of shape (n_samples, ) if method = 'int'
              ndarray of shape (n_samples, n_unique_lengths) if method = 'ohe'
              int if only sequence provided

    """

    # input handling
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
    
    # list of protein/peptide lengths
    all_len = []
    for seq in X:
        if type(seq) != str:
            seq = str(seq.seq)
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data type must be string!')
            
        all_len.append(len(seq))
    
    # one-hot-encoded lengths of proteins/peptides
    len_span = max(all_len) - min(all_len)
    len_ohe = np.zeros((len(X), len_span+1))
    len_unique = np.arange(min(all_len), max(all_len)+1)

    # fill array with ones based on sequence lengths (columns)
    for i in range(len(all_len)):    
        seq_len = all_len[i]
        j = np.where(len_unique == seq_len)[0][0]
        len_ohe[i,j] = 1

    # delete columns with all zeros
    zero_cols = np.where(~len_ohe.any(axis=0))[0]
    len_ohe = np.delete(len_ohe, zero_cols, axis=1)

    if method == 'int':
        if len(all_len) == 1:
            return all_len[0]
        else:
            return np.asarray(all_len)
    elif method == 'ohe':
        return len_ohe


def composition(X, method='relative', start=1, end=None, round_fraction=3):
    """Compute the amino acid composition of proteins or peptides.

    The frequency of each of the 20 amino acids in a protein or peptide are 
    counted and returned using IUPAC's extended protein notation of 26 letters.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    method : string, default='relative'

        'absolute' : compute absolute amino acid composition
        'relative' : compute relative amino acid composition

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

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

    # input handling
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
    
    # list of amino acids (IUPAC extended)
    amino_acids = IUPAC.ExtendedIUPACProtein().letters

    # initialize empty array with shape (n_samples, 20)
    comp_arr = np.zeros((len(X), len(amino_acids)))

    # fill array with frequency of amino acids per protein/peptide
    for i, sequence in enumerate(X):
        if type(sequence) != str:
            sequence = str(sequence.seq)
        if str.isalpha(sequence) == True:
            pass
        else:
            raise ValueError('Data type must be string!')
            
        seq = [aa for aa in sequence]
        seq = seq[start-1:end] # positional information
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
        all_len = [len(X[i]) for i in range(len(X))]
        comp_rel = [np.round((comp_arr[i]/all_len[i]), round_fraction)\
                    for i in range(comp_arr.shape[0])]
        
        return pd.DataFrame(np.asarray(comp_rel), columns=amino_acids)
        

def aaindex1(X, standardize='none', start=1, end=None):
    """Compute amino acid indices from AAIndex1.

    AAindex1 ver.9.2 (release Feb, 2017) is a set of 20 numerical values
    representing different physicochemical and biological properties of amino
    acids. Currently, it contains 566 such indices, of which 553 contain no 
    NaNs. The indices will be collected for each amino acid in the sequence, 
    then averaged across the sequence. 

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    standardize : string, default='none'

        'none' : unstandardized index matrix will be returned
        'zscore' : index matrix is standardized across columns (indices) to have
                   a mean of 0 and standard deviation of 1 (unit variance).
        'minmax' : index matrix is scaled (normalized) across columns (indices)
                   to have a range of [0, 1].

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr_index1 : Pandas DataFrame of shape (n_samples, 553-566) 

    Notes
    -----

    Columns (indices) containing NaNs will be removed. Thus, the resulting index
    matrix will have a column size of 553-566.

    The returned dataframe 'arr_index1' can easily be converted into a numpy 
    array with the command 'np.asarray(arr_index1)', if desired.
    
    Only the 20 natural amino acids are used here.

    """

    amino_acids = ['A','C','D','E','F','G','H','I','K','L',
                   'M','N','P','Q','R','S','T','V','W','Y']    
    
    # input handling
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
    
    # load AAIndex1 data
    aaind1 = pd.read_csv(PATH+'aaindex1.csv')
    _aaind1 = {}
    for aa in amino_acids:
        _aaind1[aa] = aaind1[aa].to_numpy()

    # get descriptions of all 566 indices
    desc = aaind1['Description'].values

    # initialize empty array with shape (n_samples, n_aa1_indices)
    aaind_arr = np.zeros((len(X), aaind1.shape[0]))

    # fill array with mean of indices per protein/peptide
    for i, seq in enumerate(X):
        if type(seq) != str:
            seq = str(seq.seq)
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data type must be string!')
            
        seq = seq[start-1:end] # positional information
        tmp_arr = np.zeros((aaind1.shape[0], len(seq)) )

        # fill temporary array with indices for each amino acid of the sequence
        # and compute their mean across all rows
        for j, aa in enumerate(seq):
            tmp_arr[:,j] = _aaind1[aa]

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
        

def ngram_composition(X, ngram=2, method='relative', start=1, end=None):
    """Compute n-gram peptide composition.
    
    This function computes the di-, tri-, or quadpeptide composition of
    amino acid sequences. Therefore, the function argument 'ngram' can only
    take on the values 2, 3, and 4 - otherwise, it will raise a ValueError.
    
    Parameters
    ----------
    
    X : string, fasta, or a list thereof 
       
    ngram : int, default=2
        Integer denoting the desired n-gram composition.
        
        2 : dipeptide composition
        3 : tripepitde composition
        4 : quadpeptide composition
        
    method: string, default='relative'
    
        'absolute': compute absolute ngram composition
        'relative': compute relative ngram composition

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.
        
    Returns
    -------
    
    df_ngram : Pandas DataFrame of shape (n_samples, n_unique_20^ngram)
        Depending on ngram, the returned dataframe will be of size:
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
        
    # input handling
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
    
    # get ngram combinations
    aa_combo = []
    amino_acids = IUPAC.ExtendedIUPACProtein().letters

    def combo(seq, prefix, n, k):
        "Generate every possible ngram combination"
        if k == 0: 
            aa_combo.append(prefix)
            return aa_combo  

        for i in range(n): 
            newPrefix = prefix + seq[i] 
            combo(seq, newPrefix, n, k-1) 

    combo(amino_acids, "", len(amino_acids), ngram)
    
    # create dataframe with all zeros
    arr_ngram = np.zeros((len(X), len(aa_combo)))
    df_ngram = pd.DataFrame(arr_ngram, columns=aa_combo)
    
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
    for i, sequence in enumerate(X):
        if type(sequence) != str:
            sequence = str(sequence.seq)
        if str.isalpha(sequence) == True:
            pass
        else:
            raise ValueError('Data type must be string!')

        sequence = sequence[start-1:end]
        pep_comp = n_gram(sequence)
        for j in range(len(pep_comp)):
            keys = list(pep_comp)
            df_ngram[keys[j]][i] = pep_comp[keys[j]]

    # delete zero columns
    df_ngram = df_ngram.loc[:, (df_ngram!=0).any(axis=0)]
            
    if method=='absolute':
        return df_ngram
    
    elif method=='relative':
        for i in range(df_ngram.shape[0]):
            df_ngram.iloc[i,:] = df_ngram.iloc[i,:]/sum(df_ngram.iloc[i,:])
        return df_ngram


def position_enrichment(X, position, aminoacid):
    """Compute the presence of an amino acid at a specific position.

    This function returns a binary feature vector or matrix in which ones 
    indicate the presence of the given amino acid(s) at the specified 
    position(s), and zeros indicate their absence. 

    Parameters
    ----------
    
    X : string, fasta, or a list thereof
       
    position : int or list
        Integer or list of integers denoting the position(s) in the sequence.

    aminoacid : string or list
        String or list of strings indicating the amino acid(s) of interest.
        
    Returns
    -------
    
    pos : ndarray of shape (n_samples, ) or (n_samples, n_positions)
       
    """
    
    # input handling
    # input handling
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
    
    if isinstance(position, int) and isinstance(aminoacid, str):
        pos = np.zeros((len(X),))
        for a, seq in enumerate(X):
            if type(seq) != str:
                seq = str(seq.seq)
            if str.isalpha(seq) == True:
                pass
            else:
                raise ValueError('Data type must be string!')
            for i, aa in enumerate(seq):
                if i == position-1 and aa == aminoacid:
                    pos[a] = 1
        return pos
    
    elif isinstance(position, list) and isinstance(aminoacid, list):
    
        if len(position) != len(aminoacid):
            raise ValueError("Number of positions does not match number of amino acids")

        pos = np.zeros((len(X), len(position)))

        for a, seq in enumerate(X):
            if type(seq) != str:
                seq = str(seq.seq)
            if str.isalpha(seq) == True:
                pass
            else:
                raise ValueError('Data type must be string!')
                
            for i in range(len(position)):
                if seq[position[i]-1] == aminoacid[i]:
                    pos[a, i] = 1
        return pos
    
    else:
        raise ValueError("The arguments position and aminoacid must either be integer/string or lists of integers/strings.")
