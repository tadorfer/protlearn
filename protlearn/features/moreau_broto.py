# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

# default indices of AAIndex1 (Xiao et al., 2015)
default = ['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
           'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']

def moreau_broto(X, *, d=1, properties=default, start=1, end=None): 
    """Normalized Moreau-Broto autocorrelation based on AAIndex1.

    Moreau-Broto autocorrelation descriptors are defined based on the 
    distribution of AAIndex1-based amino acid properties along the sequence. All
    indices are standardized before computing the descriptors. For the exact 
    formula, please refer to the documentation at https://protlearn.readthedocs.io/. 

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
        Array containing Moreau-Broto autocorrelation descriptors.

    References
    ----------

    Moreau & Broto (1980), Autocorrelation of a topological structure: A new 
    molecular descriptor, Nouv. J. Chim. 4, 359–360.
    
    Feng, ZP., Zhang, CT. (2000): Prediction of membrane protein types based on
    the hydrophobic index of amino acids. J Protein Chem, 19: 262–275. 
    10.1023/A:1007091128394

    Lin, Z., Pan, XM. (2001): Accurate prediction of protein secondary 
    structural content. J Protein Chem, 20: 217–220. 10.1023/A:1010967008838

    Xiao et al. (2015). protr/ProtrWeb: R package and web server for generating
    various numerical representation schemes of protein sequences. 
    Bioinformatics 31 (11), 1857-1859

    Examples
    --------

    >>> from protlearn.features import moreau_broto
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> mb = moreau_broto(seqs)
    >>> mb
    array([[ 0.20075105, -0.19610635, -0.22859247,  0.40627137, -0.19552707,
             0.16803808,  0.11262517,  0.3817956 ],
           [ 0.34693551,  0.55361789,  0.12660283,  0.37790326,  0.43905717,
             0.00559139, -0.09331898,  0.36367721]])

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

    # calculate normalized Moreau-Broto 
    arr = np.zeros((len(X), len(default)))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        for j in range(len(default)):
            p = [data[j, aadict[aa]] for aa in seq]
            ac = sum([(p[i]*p[i+d]) for i in range(len(seq)-d)])
            ats = ac/(len(seq)-d) # normalizing

            arr[i,j] = ats
            
    return arr