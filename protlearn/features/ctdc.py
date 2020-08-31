# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def ctdc(X, start=1, end=None):
    """Compute CTD composition.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 39)
    
    desc : property and group order

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
        seq = seq[start-1:end] # positional information
        tmp=[]
        for cat in categories:
            g1 = sum([seq.count(aa) for aa in group1[cat]])/len(seq)
            g2 = sum([seq.count(aa) for aa in group2[cat]])/len(seq)
            g3 = sum([seq.count(aa) for aa in group3[cat]])/len(seq)
            tmp.append([g1, g2, g3])
        arr[i,:] = [j for sublist in tmp for j in sublist]
        
    return arr, desc