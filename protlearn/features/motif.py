# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import re
import numpy as np
from ..utils.validation import check_input, check_alpha

def motif(X, pattern, start=1, end=None):
    """Returns binary vector indicating presence of custom motif.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        
    pattern : string
        Represents the sequence motif.
    
    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr :  ndarray of shape (n_samples,)

    """
    
    # input handling
    X = check_input(X)

    ## convert motif to regex pattern ##
    # replace "any"
    pattern = pattern.replace('x', '.')

    # replace "either or"
    brackets = re.findall(r'\[([A-Z]+)\]', pattern)
    if brackets:
        for rep in brackets:
            s = re.sub(r'([A-Z])(?!$)', r'\1|', rep)
            s = '(?:' + s + ')'
            pattern = re.sub(rep, s, pattern)

    # remove brackets 
    pattern = pattern.replace('[', '')
    pattern = pattern.replace(']', '')

    # replace "except"
    pattern = pattern.replace('{', '[^')
    pattern = pattern.replace('}', ']')
    
    ## compute binary vector of motif presence
    arr = np.zeros((len(X),))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        seq = seq[start-1:end] # positional information
        present = re.findall(r'{}'.format(pattern), seq)
        if present:
            arr[i] = 1
            
    return arr