# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def ctdt(X, *, start=1, end=None):
    """Composition/Transition/Distribution - Transition.

    Amino acids are categorized into 3 groups based on their physicochemical 
    properties. The properties used here include hydrophobicity, normalized van 
    der Waals volume, polarity, polarizability, charge, secondary structure, and
    solvent accessibility. For hydrophobicity, we use seven different groupings 
    based on different studies, which can all be found in AAIndex1. 

    This descriptor computes the frequency of transitions between groups. For 
    instance, if the encoded sequence is '32132223311311222222', then there are
    2 transitions from groups 1 to 2 and 2 to 1. Therefore, the descriptor for 
    this particular group transition will be 2/19. Similarly, there are 3 
    transitions from groups 2 to 3 and 3 to 2, so the descriptor for this 
    transition is 3/19. 

    As with CTDC, the dimensionality of this feature is 39.

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

    arr :  ndarray of shape (n_samples, 39)
        Array containing group transitions of physicochemical properties.

    desc : list of length 39
        Order of transition groups corresponding to columns in arr.

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

    >>> from protlearn.features import ctdt
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> t, desc = ctdt(seqs)
    >>> t
    array([[0.        , 0.        , 0.25      , 0.5       , 0.        ,
            0.        , 0.25      , 0.5       , 0.        , 0.        ,
            0.5       , 0.        , 0.25      , 0.25      , 0.        ,
            0.25      , 0.25      , 0.25      , 0.25      , 0.        ,
            0.25      , 0.        , 0.25      , 0.5       , 0.        ,
            0.25      , 0.25      , 0.        , 0.25      , 0.5       ,
            0.5       , 0.        , 0.        , 0.25      , 0.        ,
            0.        , 0.5       , 0.25      , 0.        ],
           [0.16666667, 0.33333333, 0.16666667, 0.16666667, 0.        ,
            0.        , 0.16666667, 0.        , 0.16666667, 0.16666667,
            0.        , 0.16666667, 0.5       , 0.        , 0.16666667,
            0.16666667, 0.        , 0.16666667, 0.16666667, 0.33333333,
            0.16666667, 0.16666667, 0.16666667, 0.16666667, 0.16666667,
            0.        , 0.16666667, 0.33333333, 0.        , 0.33333333,
            0.16666667, 0.16666667, 0.        , 0.        , 0.33333333,
            0.        , 0.        , 0.16666667, 0.16666667]])

    >>> desc
    ['Hydrophobicity_ARGP820101-T1221',
     'Hydrophobicity_ARGP820101-T1331',
     'Hydrophobicity_ARGP820101-T2332',
     'Hydrophobicity_CASG920101-T1221',
     'Hydrophobicity_CASG920101-T1331',
     'Hydrophobicity_CASG920101-T2332',
     'Hydrophobicity_ENGD860101-T1221',
     'Hydrophobicity_ENGD860101-T1331',
     'Hydrophobicity_ENGD860101-T2332',
     'Hydrophobicity_FASG890101-T1221',
     'Hydrophobicity_FASG890101-T1331',
     'Hydrophobicity_FASG890101-T2332',
     'Hydrophobicity_PONP930101-T1221',
     'Hydrophobicity_PONP930101-T1331',
     'Hydrophobicity_PONP930101-T2332',
     'Hydrophobicity_PRAM900101-T1221',
     'Hydrophobicity_PRAM900101-T1331',
     'Hydrophobicity_PRAM900101-T2332',
     'Hydrophobicity_ZIMJ680101-T1221',
     'Hydrophobicity_ZIMJ680101-T1331',
     'Hydrophobicity_ZIMJ680101-T2332',
     'Normalized van der Waals Volume-T1221',
     'Normalized van der Waals Volume-T1331',
     'Normalized van der Waals Volume-T2332',
     'Polarity-T1221',
     'Polarity-T1331',
     'Polarity-T2332',
     'Polarizability-T1221',
     'Polarizability-T1331',
     'Polarizability-T2332',
     'Charge-T1221',
     'Charge-T1331',
     'Charge-T2332',
     'Secondary structure-T1221',
     'Secondary structure-T1331',
     'Secondary structure-T2332',
     'Solvent accessibility-T1221',
     'Solvent accessibility-T1331',
     'Solvent accessibility-T2332']

    """

    # input handling
    X = check_input(X)

    # load data
    df = pd.read_csv(PATH+'ctd.csv')

    # get groups
    categories = list(df.Category)
    desc = [cat+'-T{}'.format(i) for cat in categories for i in ['1221','1331','2332']]
    group1 = {df['Category'][i]: df['Group1'][i] for i in range(df.shape[0])}
    group2 = {df['Category'][i]: df['Group2'][i] for i in range(df.shape[0])}
    group3 = {df['Category'][i]: df['Group3'][i] for i in range(df.shape[0])}

    # compute CTD transition
    arr = np.zeros((len(X), 39))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        cnts = []
        for cat in categories:
            # convert sequence to groups (integers)
            g1 = {aa: '1' for aa in group1[cat]}
            g2 = {aa: '2' for aa in group2[cat]}
            g3 = {aa: '3' for aa in group3[cat]}
            g = {**g1, **g2, **g3}
            conv = ''.join([g[aa] for aa in seq])
            t1221 = (conv.count('12') + conv.count('21'))/(len(seq)-1)
            t1331 = (conv.count('13') + conv.count('31'))/(len(seq)-1)
            t2332 = (conv.count('23') + conv.count('32'))/(len(seq)-1)
            cnts.append([t1221, t1331, t2332])
        arr[i,:] = [j for sublist in cnts for j in sublist]
        
    return arr, desc