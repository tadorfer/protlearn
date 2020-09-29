# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def atc(X, *, method='relative', start=1, end=None):
    """Atomic and bond composition.
    
    This function returns the sum of atomic and bond compositions for each
    amino acid sequence. The atomic features are comprised of five atoms
    (C, H, N, O, and S), and the bond features are comprised of total bonds,
    single bonds, and double bonds.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
    
    method : string, default='relative'
        'absolute': absolute atomic composition
        'relative': relative atomic composition

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr_atoms :  ndarray of shape (n_samples, 5)
        Array containing atomic compositions.
    
    arr_bonds : ndarray of shape (n_samples, 3)
        Array containing bond compositions.
    
    Notes
    -----
    
    The 'method' argument only applies to the atomic composition, not
    the bond composition.

    References
    ----------

    Kumar, R., Chaudhary, K., Singh Chauhan, J. et al. An in silico platform for
    predicting, screening and designing of antihypertensive peptides. Sci Rep 5,
    12512 (2015). https://doi.org/10.1038/srep12512

    Examples
    --------

    >>> from protlearn.features import atc
    >>> seqs = ['ARKLY', 'EERKPGL', 'AAAAAALY']
    >>> atoms, bonds = atc(seqs)
    >>> atoms
   array([[0.27522936, 0.5412844 , 0.08256881, 0.10091743, 0.],
          [0.25547445, 0.53284672, 0.08029197, 0.13138686, 0.],
          [0.26612903, 0.53225806, 0.06451613, 0.13709677, 0.]])
    >>> bonds
    array([[105.,  96.,   9.],
           [131., 121.,  10.],
           [117., 106.,  11.]])

    """

    # input handling
    X = check_input(X)

    # load data
    df = pd.read_csv(PATH+'atc.csv')
    data = np.asarray(df.iloc[:,1:])

    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    atc_dict = {amino_acids[i]: data[i,:] for i in range(20)}

    # compute atomic and bond composition
    arr_atoms = np.zeros((len(X), 5))
    arr_bonds = np.zeros((len(X), 3))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids 
        seq = seq[start-1:end] # positional information
        arr_atoms[i,:] = sum([atc_dict[aa][:5] for aa in seq]) 
        arr_bonds[i,:] = sum([atc_dict[aa][5:] for aa in seq]) 
    
    if method == 'absolute':
        return arr_atoms, arr_bonds
    
    elif method == 'relative':
        arr_atoms = [arr_atoms[i]/sum(arr_atoms[i]) for i in range(len(X))]
        return np.array(arr_atoms), arr_bonds