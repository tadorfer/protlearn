# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def distribution(seq, g):
    percentiles = [.25, .5, .75, 1]
    cnt = seq.count(g)
    if cnt:
        inds = [int(cnt*p) for p in percentiles]
        inds.insert(0, 1)
        pp = []
        for i, ind in enumerate(inds):
            pp.append([i for i, n in enumerate(seq) if n == g][ind-1]+1)
        return [i/len(seq)*100 for i in pp]
    else:
        return np.zeros((5,)).tolist()

def ctdd(X, start=1, end=None):
    """Compute CTD distribution.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 195)

    desc : order of distribution groups

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
        seq = seq[start-1:end] # positional information
        d = []
        for cat in categories:
            # convert sequence to groups (integers)
            conv = []
            for aa in seq:
                if aa in group1[cat]:
                    conv.append('1')
                elif aa in group2[cat]:
                    conv.append('2')
                elif aa in group3[cat]:
                    conv.append('3')
            conv = ''.join(conv) 
            d1 = distribution(conv, '1')
            d2 = distribution(conv, '2')
            d3 = distribution(conv, '3')
            d.append(i for j in [d1, d2, d3] for i in j)
        arr[i,:] = [j for sublist in d for j in sublist]
        
    return arr, desc