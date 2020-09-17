# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
from ..utils.validation import check_input, check_alpha, check_natural

def aac(X, method='relative', start=1, end=None):
    """Amino acid composition.

    This function returns the frequency of amino acids for each sequence in the
    dataset. 

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.

    method : string, default='relative'
        'absolute' : absolute amino acid composition
        'relative' : relative amino acid composition

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    aac :  ndarray of shape (n_samples, n_unique_amino_acids)
        Array containing the amino acid composition.
    
    amino_acids : amino acid order of aac array
        Corresponds to the columns in aac.

    Examples
    --------

    >>> from protlearn.features import aac
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> comp, aa = integer_encode(seqs)
    >>> comp
    array([[0.2       , 0.        , 0.        , 0.2       , 0.2       ,
            0.        , 0.2       , 0.2       ],
           [0.        , 0.28571429, 0.14285714, 0.14285714, 0.14285714,
            0.14285714, 0.14285714, 0.        ]])
    >>> aa
    'AEGKLPRY'

    Note that columns containing all zeros have been removed from the final 
    array.

    """

    # input handling
    X = check_input(X)
    
    # list of amino acids (IUPAC extended)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    # initialize empty array with shape (n_samples, 26)
    aac = np.zeros((len(X), len(amino_acids)))

    # compute AAC
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical
        check_natural(seq) # check for unnatural amino acids  
        seq = seq[start-1:end] # positional information
        counts = Counter(seq)
        indices = [amino_acids.index(k) for k in counts.keys()]
        aac[i,indices] = list(counts.values())

    # delete zero columns
    cols_zeros = np.where(~aac.any(axis=0))[0]
    aac = np.delete(aac, cols_zeros, axis=1)
    amino_acids = [i for j, i in enumerate(amino_acids) if j not in cols_zeros]
    amino_acids = ''.join(amino_acids)

    if method == 'absolute':
        return aac, amino_acids

    elif method == 'relative':
        l = [len(seq) for seq in X]
        aac_rel = [(aac[i]/l[i]) for i in range(aac.shape[0])]
        return np.stack(aac_rel), amino_acids