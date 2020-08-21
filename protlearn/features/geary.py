# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from Bio.Alphabet import IUPAC
from utils.validation import check_input
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

# default indices of AAIndex1 (Xiao et al., 2015)
default = ['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
           'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']

def geary(X, d=1, properties=default, start=1, end=None): 
    """Compute Geary's C based on AAIndex1 features.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
    
    properties : list
        List of strings denoting AAIndex1 indices.
        
    d : int, default=1
        Represents the lag. Must be smaller than sequence length.
        Maximum: 30.
    
    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, n_properties)

    """

    # input handling
    X = check_input(X)
    min_len = min([len(seq) for seq in X])
    if d >= 30:
        raise ValueError('Maximum lag parameter is 30!')
    if d >= min_len:
        raise ValueError('Lag parameter d must be smaller than sequence length!')
    
    # load data
    df = pd.read_csv(PATH+'aaindex1.csv').set_index('Description')
    df = df.reindex(sorted(df.columns), axis=1)
    data = np.asarray(df.loc[default])

    # list of amino acids (IUPAC standard)
    amino_acids = IUPAC.IUPACProtein.letters
    aadict = {amino_acids[i]: i for i in range(20)}

    # standardization
    for i in range(data.shape[0]):
        data[i,:] = [(j-np.mean(data[i,:]))/np.std(data[i,:]) for j in data[i,:]]

    # calculate Geary's C
    arr = np.zeros((len(X), len(default)))
    for i, seq in enumerate(X):
        # check that input is alphabetical
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')

        seq = seq[start-1:end] # positional information
        eq1 = 1/(2*(len(seq)-d))
        for j in range(len(default)):
            p = [data[j, aadict[aa]] for aa in seq]
            p_prime = sum(p)/len(seq)
            eq2 = sum([(p[i]-p[i+d])**2 for i in range(len(seq)-d)])
            eq3 = sum([(p[i]-p_prime)**2 for i in range(len(seq))])

            arr[i,j] = (eq1*eq2)/((1/(len(seq)-1))*eq3)
            
    return arr