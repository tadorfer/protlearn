# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def ctdt(X, start=1, end=None):
    """Compute CTD transition.

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

    desc : order of transition groups

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
        seq = seq[start-1:end] # positional information
        cnts = []
        for cat in categories:
            # convert sequence to groups (integers)
            tmp = []
            for aa in seq:
                if aa in group1[cat]:
                    tmp.append('1')
                elif aa in group2[cat]:
                    tmp.append('2')
                elif aa in group3[cat]:
                    tmp.append('3')
            tmp = ''.join(tmp)
            t1221 = (tmp.count('12') + tmp.count('21'))/(len(seq)-1)
            t1331 = (tmp.count('13') + tmp.count('31'))/(len(seq)-1)
            t2332 = (tmp.count('23') + tmp.count('32'))/(len(seq)-1)
            cnts.append([t1221, t1331, t2332])
        arr[i,:] = [j for sublist in cnts for j in sublist]
        
    return arr, desc