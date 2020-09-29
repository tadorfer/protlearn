# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

# default indices of AAIndex1 (Xiao et al., 2015)
default = ['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
           'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']

def geary(X, *, d=1, properties=default, start=1, end=None): 
    """Geary's C based on AAIndex1.

    Geary's C autocorrelation descriptors are defined based on the distribution 
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
        Array containing Geary's C autocorrelation descriptors.

    References
    ----------

    Geary, R. (1954). The Contiguity Ratio and Statistical Mapping. The 
    Incorporated Statistician, 5(3), 115-146. doi:10.2307/2986645

    Jeffers, J. (1973). A Basic Subroutine for Geary's Contiguity Ratio. Journal
    of the Royal Statistical Society. Series D (The Statistician), 22(4), 
    299-302. doi:10.2307/2986827

    Sokal, RR., Thomson, BA. (2006). Population structure inferred by local spatial 
    autocorrelation: an example from an Amerindian tribal population. Am J Phys
    Anthropol, 129: 121–131. 10.1002/ajpa.20250

    Xiao et al. (2015). protr/ProtrWeb: R package and web server for generating
    various numerical representation schemes of protein sequences. 
    Bioinformatics 31 (11), 1857-1859

    Examples
    --------

    >>> from protlearn.features import geary
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> gearyC = geary(seqs)
    >>> gearyC
    array([[0.52746275, 1.12898944, 0.94222955, 0.39077186, 0.96444569,
            0.66346012, 0.87481962, 0.32546227],
           [0.65656058, 0.95397893, 0.87962853, 0.70972353, 0.65407555,
            0.96823847, 1.01949384, 0.21073089]])

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

    # calculate Geary's C
    arr = np.zeros((len(X), len(default)))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        eq1 = 1/(2*(len(seq)-d))
        for j in range(len(default)):
            p = [data[j, aadict[aa]] for aa in seq]
            p_prime = sum(p)/len(seq)
            eq2 = sum([(p[i]-p[i+d])**2 for i in range(len(seq)-d)])
            eq3 = sum([(p[i]-p_prime)**2 for i in range(len(seq))])

            arr[i,j] = (eq1*eq2)/((1/(len(seq)-1))*eq3)
            
    return arr