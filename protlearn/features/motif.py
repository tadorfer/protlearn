# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import re
import numpy as np
from ..utils.validation import check_input, check_alpha

def motif(X, pattern, *, start=1, end=None):
    """Sequence motifs.

    This function returns a binary vector indicating the presence of a specified 
    amino acid sequence motif.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
        
    pattern : string
        Represents the sequence motif.
        x --> any amino acid
        [XY] --> X or Y
        {X} --> any amino acid except X
    
    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr :  ndarray of shape (n_samples,)
        Binary vector indicating the presence of the motif in sequences.   

    Examples
    --------

    >>> from protlearn.features import motif
    >>> seqs = ['AARKYLL', 'LELCDPGPG', 'RAAANCDD']  
    >>> pattern1 = pattern = 'AAx[KC]'
    >>> m1 = motif(seqs, pattern1)
    >>> m1
    array([1., 0., 1.])
    >>> pattern2 = 'xxC[DA]xx{Y}'
    >>> m2 = motif(seqs, pattern2)
    >>> m2
    array([0., 1., 0.]) 

    Based on the example above, 'pattern1' is interpreted as follows:
    Two consecutive amino acids 'A', followed by any amino acid, followed by
    either a 'K' or a 'C'. 

    Likewise, pattern2 is interpreted as follows:
    Any two consecutive amino acids, followed by a 'C', followed by either a 'D'
    or an 'A', followed by any two amino acids, followed by any amino acid
    except 'Y'.

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