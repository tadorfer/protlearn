# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

# default indices of AAIndex1 (Xiao et al., 2015)
default = ['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
           'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']

def moran(X, *, d=1, properties=default, start=1, end=None): 
    """Moran's I based on AAIndex1.

    Moran's I autocorrelation descriptors are defined based on the distribution 
    of AAIndex1-based amino acid properties along the sequence. All indices are 
    standardized before computing the descriptors. For the exact formula, please
    refer to the documentation at https://protlearn.readthedocs.io/. 

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
    
    properties : list
        List of strings denoting AAIndex1 indices.
        
    d : int, default=1
        Represents the lag. Must be smaller than sequence length. Maximum: 30.
    
    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr : ndarray of shape (n_samples, n_properties)   
        Array containing Moran's I autocorrelation descriptors.

    References
    ----------

    Moran, P. (1950). Notes on Continuous Stochastic Phenomena. Biometrika, 
    37(1/2), 17-23. doi:10.2307/2332142
    
    Horne, DS. (1988): Prediction of protein helix content from an 
    autocorrelation analysis of sequence hydrophobicities. Biopolymers 27, 
    451–477. 10.1002/bip.360270308

    Li et al. (2007). Beyond Moran's I: testing for spatial
    dependence based on the spatial autoregressive model. Geogr. Anal. 39, 
    357–375.

    Xiao et al. (2015). protr/ProtrWeb: R package and web server for generating
    various numerical representation schemes of protein sequences. 
    Bioinformatics 31 (11), 1857-1859

    Examples
    --------

    >>> from protlearn.features import moran
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> moranI = moran(seqs)
    >>> moranI
    array([[ 0.40090094, -0.31240708, -0.44083728,  0.26720303, -0.45198768,
            -0.14684112, -0.05212843,  0.33703981],
           [-0.0588976 , -0.36033526,  0.13170834,  0.18317369,  0.3884609 ,
            -0.00724234, -0.19231646,  0.61711506]])

    """

    # input handling
    X = check_input(X)
    min_len = min([len(seq) for seq in X])
    if d > 30:
        raise ValueError('Maximum lag parameter is 30!')
    if d >= min_len:
        raise ValueError('Lag parameter d must be smaller than sequence length!')
    
    # load data
    df = pd.read_csv(PATH+'aaindex1.csv').set_index('Description')
    df = df.reindex(sorted(df.columns), axis=1)
    data = np.asarray(df.loc[default])

    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aadict = {amino_acids[i]: i for i in range(20)}

    # standardization
    for i in range(data.shape[0]):
        data[i,:] = [(j-np.mean(data[i,:]))/np.std(data[i,:]) for j in data[i,:]]

    # calculate Moran's I
    arr = np.zeros((len(X), len(default)))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        eq1 = 1/(len(seq)-d)  
        for j in range(len(default)):
            p = [data[j, aadict[aa]] for aa in seq]
            p_prime = sum(p)/len(seq)
            eq2 = sum([(p[i]-p_prime)*(p[i+d]-p_prime) for i in range(len(seq)-d)])
            eq3 = sum([(p[i]-p_prime)**2 for i in range(len(seq))])

            arr[i,j] = (eq1*eq2)/((1/len(seq))*eq3)
            
    return arr