# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def distribution(seq, g):
    percentiles = [.25, .5, .75, 1]
    cnt = seq.count(g)
    if cnt:
        inds = [int(cnt*p) for p in percentiles]
        inds.insert(0, 1)
        pp = []
        for ind in inds:
            pp.append([i for i, n in enumerate(seq) if n == g][ind-1]+1)
        return [i/len(seq)*100 for i in pp]
    else:
        return np.zeros((5,)).tolist()

def ctdd(X, *, start=1, end=None):
    """Composition/Transition/Distribution - Distribution.

    Amino acids are categorized into 3 groups based on their physicochemical 
    properties. The properties used here include hydrophobicity, normalized van 
    der Waals volume, polarity, polarizability, charge, secondary structure, and
    solvent accessibility. For hydrophobicity, we use seven different groupings 
    based on different studies, which can all be found in AAIndex1. 
    
    There are five distribution descriptors for each physicochemical property
    and they are the position percents in the whole sequence for the first 
    residue, 25% residues, 50% residues, 75% residues, and 100% residues for a 
    certain encoded class. For instance, if the encoded sequence is 
    '32132223311311222222', then there are 10 residues encoded as 2. The 
    positions for the first residue 2, the 2nd residue 2 (25% * 10 = 2), the 5th
    2 residue (50% * 10 = 5), the 7th 2 (75% * 10 = 7) and the 10th residue 2 
    (100% * 10) in the encoded sequence are 2, 5, 15, 17, 20, so that the 
    distribution descriptors for 2 are: 10.0 (2/20 * 100), 25.0 (5/20 * 100), 
    75.0 (15/20 * 100), 85.0 (17/20 * 100), 100.0 (20/20 * 100).

    Since there are 13 physicochemical properties, 3 groups (see ctdc), and 5 
    distribution descriptors, the total dimensionality of this feature is 195.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 195)
        Array containing grouped distribution of physicochemical properties.

    desc : list of length 195
        Order of distribution groups corresponding to columns in arr.

    References
    ----------

    Dubchak, I., Muchnik, I., Holbrook, S. R. & Kim, S.-H. Prediction of protein
    folding class using global description of amino acid sequence. Proceedings 
    of the National Academy of Sciences 92, 8700–8704 (1995).

    Dubchak, I., Muchnik, I., Mayor, C., Dralyuk, I. & Kim, S.-H. Recognition of
    a protein fold in the context of the scop classification. Proteins: 
    Structure, Function, and Bioinformatics 35, 401–407 (1999).

    Examples
    --------

    >>> from protlearn.features import ctdd
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> d, desc = ctdd(seqs)
    >>> d.shape
    (2, 195)
    >>> len(desc)
    195

    """

    # input handling
    X = check_input(X)

    # load data
    df = pd.read_csv(PATH+'ctd.csv')

    # get groups
    categories = list(df.Category)
    groups = ['1', '2', '3']
    percentiles = ['0', '25', '50', '75', '100']
    desc = [cat+'-G{}'.format(i) for cat in categories for i in groups]
    desc = [d+'D{}'.format(i) for d in desc for i in percentiles]
    group1 = {df['Category'][i]: df['Group1'][i] for i in range(df.shape[0])}
    group2 = {df['Category'][i]: df['Group2'][i] for i in range(df.shape[0])}
    group3 = {df['Category'][i]: df['Group3'][i] for i in range(df.shape[0])}

    # compute CTD distribution
    arr = np.zeros((len(X), 195))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        d = []
        for cat in categories:
            # convert sequence to groups (integers)
            g1 = {aa: '1' for aa in group1[cat]}
            g2 = {aa: '2' for aa in group2[cat]}
            g3 = {aa: '3' for aa in group3[cat]}
            g = {**g1, **g2, **g3}
            conv = ''.join([g[aa] for aa in seq])
            d1 = distribution(conv, '1')
            d2 = distribution(conv, '2')
            d3 = distribution(conv, '3')
            d.append(i for j in [d1, d2, d3] for i in j)
        arr[i,:] = [j for sublist in d for j in sublist]
        
    return arr, desc