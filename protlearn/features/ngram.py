# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
import numpy as np
import pandas as pd
from collections import Counter
from Bio.Alphabet import IUPAC
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pkg_resources

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