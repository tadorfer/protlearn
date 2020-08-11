# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
from Bio.Alphabet import IUPAC
from utils.validation import check_input

def aac(X, method='relative', start=1, end=None):
    """Compute the amino acid composition of proteins or peptides.

    The frequency of each of amino acids in a protein or peptide are counted 
    and returned using IUPAC's extended protein notation of 26 letters.

    Parameters
    ----------

    X : string, fasta, or a list thereof 

    method : string, default='relative'

        'absolute' : compute absolute amino acid composition
        'relative' : compute relative amino acid composition

    start : int, default=1
        Determines the starting point of the amino acid sequence.

    end : int, default=None
        Determines the end point of the amino acid sequence.

    Returns
    -------

    aac :  ndarray of shape (n_samples, n_unique_amino_acids)
    
    amino_acids : amino acid order of aac array

    Notes
    -----

    Start and end positions can be specified to determine the composition of
    a particular amino acid stretch within the protein.

    """

    # input handling
    X = check_input(X)
    
    # list of amino acids (IUPAC extended)
    amino_acids = IUPAC.ExtendedIUPACProtein().letters

    # initialize empty array with shape (n_samples, 26)
    aac = np.zeros((len(X), len(amino_acids)))

    # compute AAC
    for i, seq in enumerate(X):
        # if fasta, get sequence as string
        if type(seq) != str:
            seq = str(seq.seq)
            
        # check that input is alphabeticalb
        if str.isalpha(seq) == True:
            pass
        else:
            raise TypeError('Data type must be string!')
            
        seq = seq[start-1:end] # positional information
        counts = Counter(seq)
        indices = [amino_acids.index(k) for k in counts.keys()]
        aac[i,indices] = list(counts.values())

    # delete zero columns
    cols_zeros = np.where(~aac.any(axis=0))[0]
    aac = np.delete(aac, cols_zeros, axis=1)
    amino_acids = [i for j, i in enumerate(amino_acids) if j not in cols_zeros]

    if method == 'absolute':
        return aac, amino_acids

    elif method == 'relative':
        l = [len(seq) for seq in X]
        aac_rel = [(aac[i]/l[i]) for i in range(aac.shape[0])]
        return np.stack(aac_rel), amino_acids