# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from ..utils.validation import check_input, check_alpha

def remove_unnatural(X):
    """Remove sequences containing unnatural amino acids.

    This function removes sequences containing amino acids other than the 20 
    natural ones.

    Parameters
    ----------

    X : string, fasta, or a list thereof
        Dataset of amino acid sequences.

    Returns
    -------

    Y : list of length n_samples minus the number of sequences containing 
        unnatural amino acids
        Dataset containing only sequences comprised of natural amino acids.

    Examples
    --------

    >>> from protlearn.preprocessing import remove_unnatural
    >>> seqs = ['ARKLY', 'EERNPJAB', 'QEPGPGLLLK']
    >>> seqs = remove_unnatural(seqs)
    >>> seqs
    ['ARKLY', 'QEPGPGLLLK']

    """
    
    # input handling 
    X = check_input(X)

    # remove sequences with unnatural amino acids
    amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    indices = []
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical      
        # get indices of sequences with unnatural amino acids
        for aa in seq:
            if aa in amino_acids:
                pass
            else:
                indices.append(i)
                break

    # remove sequences by indices
    Y = [i for j, i in enumerate(X) if j not in set(indices)]
    
    return Y
