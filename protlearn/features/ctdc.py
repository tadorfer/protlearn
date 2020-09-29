# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def ctdc(X, *, start=1, end=None):
    """Composition/Transition/Distribution - Composition.

    Amino acids are categorized into 3 groups based on their physicochemical 
    properties. The properties used here include hydrophobicity, normalized van 
    der Waals volume, polarity, polarizability, charge, secondary structure, and
    solvent accessibility. For hydrophobicity, we use seven different groupings 
    based on different studies, which can all be found in AAIndex1. 
    
    After grouping, the frequency of each class will be calculated for each 
    physicochemical property per sequence in the dataset. For instance, the 
    sequence 'ARKLY' translates to '23311' with respect to polarity groups. 
    Thus, for polarity, the outcome will be P1 = 2/5, P2 = 1/5, and P3 = 2/5. 
    
    As there are 13 different properties and 3 groups for each, the total 
    dimension of this descriptor is 39. 

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

    arr : ndarray of shape (n_samples, 39)
        Array containing grouped composition of physicochemical properties.
    
    desc : list of length 39
        Descriptor properties corresponding to columns in arr.

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

    >>> from protlearn.features import ctdc
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> c, desc = ctdc(seqs)
    >>> c
    array([[0.        , 0.6       , 0.4       , 0.4       , 0.6       ,
            0.        , 0.6       , 0.2       , 0.2       , 0.4       ,
            0.        , 0.6       , 0.2       , 0.4       , 0.4       ,
            0.4       , 0.4       , 0.2       , 0.4       , 0.2       ,
            0.4       , 0.2       , 0.2       , 0.6       , 0.4       ,
            0.2       , 0.4       , 0.2       , 0.2       , 0.6       ,
            0.4       , 0.6       , 0.        , 0.8       , 0.2       ,
            0.        , 0.4       , 0.4       , 0.2       ],
           [0.42857143, 0.28571429, 0.28571429, 0.85714286, 0.14285714,
            0.        , 0.71428571, 0.14285714, 0.14285714, 0.57142857,
            0.28571429, 0.14285714, 0.57142857, 0.28571429, 0.14285714,
            0.57142857, 0.28571429, 0.14285714, 0.57142857, 0.14285714,
            0.28571429, 0.28571429, 0.42857143, 0.28571429, 0.14285714,
            0.28571429, 0.57142857, 0.14285714, 0.57142857, 0.28571429,
            0.28571429, 0.42857143, 0.28571429, 0.71428571, 0.        ,
            0.28571429, 0.28571429, 0.57142857, 0.14285714]])
    >>> desc
    ['Hydrophobicity_ARGP820101-G1',
     'Hydrophobicity_ARGP820101-G2',
     'Hydrophobicity_ARGP820101-G3',
     'Hydrophobicity_CASG920101-G1',
     'Hydrophobicity_CASG920101-G2',
     'Hydrophobicity_CASG920101-G3',
     'Hydrophobicity_ENGD860101-G1',
     'Hydrophobicity_ENGD860101-G2',
     'Hydrophobicity_ENGD860101-G3',
     'Hydrophobicity_FASG890101-G1',
     'Hydrophobicity_FASG890101-G2',
     'Hydrophobicity_FASG890101-G3',
     'Hydrophobicity_PONP930101-G1',
     'Hydrophobicity_PONP930101-G2',
     'Hydrophobicity_PONP930101-G3',
     'Hydrophobicity_PRAM900101-G1',
     'Hydrophobicity_PRAM900101-G2',
     'Hydrophobicity_PRAM900101-G3',
     'Hydrophobicity_ZIMJ680101-G1',
     'Hydrophobicity_ZIMJ680101-G2',
     'Hydrophobicity_ZIMJ680101-G3',
     'Normalized van der Waals Volume-G1',
     'Normalized van der Waals Volume-G2',
     'Normalized van der Waals Volume-G3',
     'Polarity-G1',
     'Polarity-G2',
     'Polarity-G3',
     'Polarizability-G1',
     'Polarizability-G2',
     'Polarizability-G3',
     'Charge-G1',
     'Charge-G2',
     'Charge-G3',
     'Secondary structure-G1',
     'Secondary structure-G2',
     'Secondary structure-G3',
     'Solvent accessibility-G1',
     'Solvent accessibility-G2',
     'Solvent accessibility-G3']

    """

    # input handling
    X = check_input(X)

    # load data
    df = pd.read_csv(PATH+'ctd.csv')

    # get groups
    categories = list(df.Category)
    desc = [cat+'-G{}'.format(i) for cat in categories for i in range(1,4)]
    group1 = {df['Category'][i]: df['Group1'][i] for i in range(df.shape[0])}
    group2 = {df['Category'][i]: df['Group2'][i] for i in range(df.shape[0])}
    group3 = {df['Category'][i]: df['Group3'][i] for i in range(df.shape[0])}

    # compute CTD composition
    arr = np.zeros((len(X), 39))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        tmp=[]
        for cat in categories:
            g1 = sum([seq.count(aa) for aa in group1[cat]])/len(seq)
            g2 = sum([seq.count(aa) for aa in group2[cat]])/len(seq)
            g3 = sum([seq.count(aa) for aa in group3[cat]])/len(seq)
            tmp.append([g1, g2, g3])
        arr[i,:] = [j for sublist in tmp for j in sublist]
        
    return arr, desc