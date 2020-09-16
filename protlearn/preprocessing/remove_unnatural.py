# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from ..utils.validation import check_input, check_alpha

def remove_unnatural(X):
    """Remove sequences containing unnatural amino acids.

    Parameters
    ----------

    X : string, fasta, or a list thereof

    Returns
    -----

    Y : list of length n_samples minus the number of sequences containing 
        unnatural amino acids

    """
    
    # input handling 
    X = check_input(X)

    # remove sequences with unnatural amino acids
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    indices = []
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical      
        for aa in seq:
            if aa in set(amino_acids):
                pass
            else:
                indices.append(i)

    Y = [i for j, i in enumerate(X) if j not in set(indices)]
    
    return Y