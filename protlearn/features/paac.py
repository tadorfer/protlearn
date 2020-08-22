# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from Bio.Alphabet import IUPAC
from ..utils.validation import check_input
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def paac(X, lambda_=1, w=.05, start=1, end=None):
    """Compute pseudo amino acid composition.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
    
    lambda_ : int, default=1
        Counted rank (tier) of the correlation along an amino acid sequence.
        
    w : float, default=.05
        Weight factor.

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 20+lambda_)
    
    amino_acids : amino acid order of paac array

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
        desc.append('lambda' + str(n))
    
    # normalization
    for i in range(data.shape[1]):
        mean = np.mean(data[:,i])
        denom = np.sqrt(sum([(j-mean)**2 for j in data[:,i]])/20)
        data[:,i] = [(j-mean)/denom for j in data[:,i]]

    aa_dict = {amino_acids[i]: data[i,:] for i in range(20)}

    # correlation function
    def corr(aa1, aa2, data):
        return sum([(aa_dict[aa1][i] - aa_dict[aa2][i])**2 \
                    for i in range(data.shape[1])]) / data.shape[1]

    # computing pseudo amino acid composition
    arr = np.zeros((len(X),len(amino_acids)+lambda_))
    for i, seq in enumerate(X):
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')

        seq = seq[start-1:end] # positional information
        theta = []
        for n in range(1, lambda_+1):
            theta.append(
                sum([corr(seq[j], seq[j+n], data) \
                     for j in range(len(seq)-n)])/(len(seq)-n))
        
        cnt = Counter(seq)
        pse_aac = [cnt[aa] / (1+w * sum(theta)) for aa in amino_acids]
        pse_aac = pse_aac + [(w*j) / (1+w * sum(theta)) for j in theta]
        arr[i,:] = pse_aac

    # delete zero columns
    cols_zeros = np.where(~arr.any(axis=0))[0]
    arr = np.delete(arr, cols_zeros, axis=1)
    desc = [i for j, i in enumerate(desc) if j not in cols_zeros]

    return arr, desc