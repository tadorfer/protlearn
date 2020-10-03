# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from ..utils.validation import check_input, check_alpha, check_natural

def binary(X, *, padding=True, start=1, end=None):
    """Binary profile pattern.

    This function returns the binary profile pattern for each amino acid 
    sequence in the dataset. The output array is therefore very sparse, as each
    amino acid is represented by a 20-dimensional vector with only one non-zero
    value.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
    
    padding : bool, default=True
        Pad sequences of unequal lengths with zeros at the posterior end.

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 20*seq_length)
        Array containing binary profile pattern.
    
    Notes
    -----
    
    This function is intended for proteins or peptides with equal lengths only.

    References
    ----------

    Ansari, H. R., & Raghava, G. P. (2010). Identification of conformational 
    B-cell Epitopes in an antigen from its primary sequence. Immunome research,
    6, 6. https://doi.org/10.1186/1745-7580-6-6

    Chen Z, Chen Y-Z, Wang X-F, Wang C, Yan R-X, Zhang Z (2011) Prediction of 
    Ubiquitination Sites by Using the Composition of k-Spaced Amino Acid Pairs.
    PLoS ONE 6(7): e22930. https://doi.org/10.1371/journal.pone.0022930

    Examples
    --------

    >>> from protlearn.features import binary
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> bpp = binary(seqs, padding=True)
    >>> bpp
    array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])

    """

    # input handling
    X = check_input(X)

    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_dict = {aa: i for i, aa in enumerate(amino_acids)}

    # define maximum length 
    l = [len(seq) for seq in X]
    max_len = max(l)
    if padding == False and len(l) != l.count(max_len):
        raise ValueError('Sequences must be of equal length or padded!')
        
    # compute binary profile pattern
    arr = np.zeros((len(X), len(amino_acids)*max_len))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        binary = [aa_dict[aa]+x*20 for x, aa in enumerate(seq)]
        arr[i,binary] = 1

    return arr