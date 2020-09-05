# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from ..utils.validation import check_input, check_alpha
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def atc(X, method='relative', start=1, end=None):
    """Compute atomic and bond composition.
    
    This function returns the sum of atomic and bond compositions per
    amino acid sequence. The atomic features are comprised of five atoms
    (C, H, N, O, and S), and the bond features are comprised of total bonds,
    single bonds, and double bonds.

    Parameters
    ----------

    X : string, fasta, or a list thereof 
    
    method: string, default='relative'
    
        'absolute': compute absolute atomic composition
        'relative': compute relative atomic composition

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    arr_atoms :  ndarray of shape (n_samples, 5)
    
    arr_bonds : ndarray of shape (n_samples, 3)
    
    Notes
    -----
    
    The 'method' argument only applies to the atomic composition, not
    the bond composition.

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
        seq = seq[start-1:end] # positional information
        arr_atoms[i,:] = sum([atc_dict[aa][:5] for aa in seq]) 
        arr_bonds[i,:] = sum([atc_dict[aa][5:] for aa in seq]) 
    
    if method == 'absolute':
        return arr_atoms, arr_bonds
    
    elif method == 'relative':
        arr_atoms = [arr_atoms[i]/sum(arr_atoms[i]) for i in range(len(X))]
        return arr_atoms, arr_bonds