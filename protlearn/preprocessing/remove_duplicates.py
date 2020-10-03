# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from collections import Counter
from ..utils.validation import check_input

def remove_duplicates(X, *, verbose=1):
    """Remove duplicate sequences.

    This function detects and removes duplicate sequences from the dataset.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
    
    verbose : int, default=1
        0 : no information on duplicates is printed
        1 : prints number of duplicates removed
        2 : prints duplicate sequences and number of times present

    Returns
    -------

    Y : list of length n_samples minus the number of duplicates
        Dataset containing only unique sequences.

    Examples
    --------

    >>> from protlearn.preprocessing import remove_duplicates
    >>> seqs = ['ARKLY', 'EERNPAA', 'ARKLY', 'QEPGPGLLLK']
    >>> seqs = remove_duplicates(seqs)
    >>> seqs
    ['EERNPAA', 'QEPGPGLLLK', 'ARKLY']

    """

    # input handling
    X = check_input(X)

    # remove duplicates
    Y = list(set(X))
    diff = Counter(X)-Counter(Y)
    
    # handle verbosity
    if diff:
        if verbose == 1:
            if len(diff) == 1:
                print("1 duplicate has been removed.")
            else:
                print("{} duplicates have been removed.".format(len(diff)))
        elif verbose == 2:
            for entry in diff:
                print("Removed duplicates [sequence: number of times present]:\n")
                print("{0}: {1}".format(entry, diff[entry]+1))
    else:
        print("Data contains no duplicates.")
    
    return Y