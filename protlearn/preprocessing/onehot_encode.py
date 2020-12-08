# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from ..utils.validation import check_input, check_alpha, check_natural

def onehot_encode(X):
    """One-hot encoding.
    
    This function converts amino acid sequences into their corresponding
    one-hot encoded representations. Sequences will be padded with zeros to the
    maximum sequence length so that the final output has the shape of 
    (n_samples, maximum_length, 20), where 20 is the number of natural amino 
    acids.
    
    Parameters
    ----------

    X : string, fasta, or a list thereof
        Dataset of amino acid sequences.

    Returns
    -------

    enc : ndarray of shape (n_samples, max_len, 20) 
        Contains the one-hot-encoded amino acid sequences.

    Examples
    --------

    >>> from protlearn.preprocessing import onehot_encode
    >>> seqs = ['ARKLY', 'EERNPAA', 'QEPGPGLLLK']
    >>> enc = onehot_encode(seqs, padding=True)
    >>> enc.shape
    (3, 10, 20)
    
    """

    # input handling
    X = check_input(X)
    
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_dict = {aa: i for i, aa in enumerate(amino_acids)}

    # get maximum sequence length in dataset
    max_len = max([len(seq) for seq in X])

    arr = np.zeros((len(X), max_len, len(amino_acids)))

    # one-hot encoding
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical      
        check_natural(seq) # check for unnatural amino acids
        for j, aa in enumerate(seq):
            k = aa_dict[aa]
            arr[i,j,k] = 1

    return arr