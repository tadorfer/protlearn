# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from Bio.Alphabet import IUPAC
from ..utils.validation import check_input
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def apaac(X, lambda_=30, w=.05, start=1, end=None):
    """Compute amphiphilic pseudo amino acid composition.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
    
    lambda_ : int, default=30
        Counted rank (tier) of the correlation along an amino acid sequence.
        
    w : float, default=.05
        Weight factor.

    start : int, default=30
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 20+2*lambda_)
    
    desc : amino acid + lambda order of apaac array

    """
    
    # input handling
    X = check_input(X)
    
    # load data
    df = pd.read_csv(PATH+'paac.csv')
    data = np.asarray(df.iloc[:,1:])

    # list of amino acids (IUPAC standard)
    amino_acids = IUPAC.IUPACProtein.letters
    desc = [aa for aa in amino_acids]

    for n in range(1, lambda_+1):
        desc.append('lambda_hphob' + str(n))
        desc.append('lambda_hphil' + str(n))
    
    # normalization
    for i in range(data.shape[1]):
        mean = np.mean(data[:,i])
        denom = np.sqrt(sum([(j-mean)**2 for j in data[:,i]])/20)
        data[:,i] = [(j-mean)/denom for j in data[:,i]]

    aa_dict = {amino_acids[i]: data[i,:] for i in range(20)}

    # computing pseudo amino acid composition
    arr = np.zeros((len(X),len(amino_acids)+2*lambda_))
    for i, seq in enumerate(X):
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')

        seq = seq[start-1:end] # positional information
        tau = []
        for n in range(1, lambda_+1):
            for j in range(2):
                tau.append(
                    sum([aa_dict[seq[k]][j]*aa_dict[seq[k+n]][j] \
                         for k in range(len(seq)-n)])/(len(seq)-n))
        
        cnt = Counter(seq)
        apse_aac = [cnt[aa] / (1+w * sum(tau)) for aa in amino_acids]
        apse_aac = apse_aac + [(w*j) / (1+w * sum(tau)) for j in tau]
        arr[i,:] = apse_aac

    # delete zero columns
    cols_zeros = np.where(~arr.any(axis=0))[0]
    arr = np.delete(arr, cols_zeros, axis=1)
    desc = [i for j, i in enumerate(desc) if j not in cols_zeros]

    return arr, desc