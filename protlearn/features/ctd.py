# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
from itertools import product
from ..utils.validation import check_input, check_alpha, check_natural

def ctd(X, *, start=1, end=None):
    """Conjoint triad descriptors.

    These descriptors were initially developed to model protein-protein
    interactions. Amino acids can be grouped into 7 different classes based on
    their dipoles and side chain volumes, which reflect their electrostatic and
    hydrophobic interactions. After grouping, these class triads are counted and
    normalized (for more information, see https://protlearn.readthedocs.io/).

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

    arr :  ndarray of shape (n_samples, 343)
        Array containing conjoint triad descriptors.

    ctd_list : list of length 343
        Unique class triads corresponding to columns in arr.

    References
    ----------

    Shen J, Zhang J, Luo X, Zhu W, Yu K, Chen K, Li Y, Jiang H (2007) Predicting
    protein-protein interactions based only on sequences information. Proc Natl
    Acad Sci USA 104: 4337 – 4341

    Examples
    --------

    >>> from protlearn.features import ctd
    >>> seqs = ['ARKKLYLYL', 'EEEERKPGL']
    >>> ctd_arr, ctd_desc = ctd(seqs)
    >>> ctd_arr.shape
    (2, 343)
    >>> len(ctd_desc)
    343

    """

    # input handling
    X = check_input(X)

    # define classes
    classes = {'A': 1, 'G': 1, 'V': 1,
               'I': 2, 'L': 2, 'F': 2, 'P': 2,
               'Y': 3, 'M': 3, 'T': 3, 'S': 3,
               'H': 4, 'N': 4, 'Q': 4, 'W': 4,
               'R': 5, 'K': 5,
               'D': 6, 'E': 6,
               'C': 7}

    # compute CTD
    ctd_list = [''.join(i) for i in product('1234567', repeat=3)]
    ctd = {ctd_list[i]: [] for i in range(len(ctd_list))}
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        seq = ''.join([str(classes[aa]) for aa in seq])
        keys = [seq[x:x+3] for x in range(len(seq)-2)]
        unq = set(keys)
        other = ctd.keys()-unq
        ctd_vals = sorted(Counter(keys).items())
        vals = [i[1] for i in ctd_vals]

        for num, j in enumerate(unq):
            if j in ctd:
                ctd[j].append(vals[num])
            
        for k in other:
            ctd[k].append(0)
    
    arr = np.array(list(ctd.values()), dtype=float).T
    
    return arr, ctd_list