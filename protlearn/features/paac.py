# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def paac(X, *, lambda_=30, w=.05, remove_zero_cols=False, start=1, end=None):
    """Pseudo amino acid composition.

    Similar to the vanilla amino acid composition, this feature characterizes 
    the protein mainly using a matrix of amino-acid frequencies, which helps 
    with dealing with proteins without significant sequence homology to other 
    proteins. However, additional information are also included in the 
    matrix to represent some local features, such as correlation between 
    residues of a certain distance. For the exact method of calculating the 
    pseudo amino acid composition, please refer to the documentation at
    https://protlearn.readthedocs.io/. 

    Parameters
    ----------

    X : string, fasta, or a list thereof 
        Dataset of amino acid sequences.
    
    lambda_ : int, default=30
        Counted rank (tier) of the correlation along an amino acid sequence.
        This parameter has to be smaller than the shortest sequence in the 
        dataset.
        
    w : float, default=.05
        Weighting factor for the sequence-order effect.

    remove_zero_cols : bool, default=False
        If true, columns containing only zeros will be deleted.

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr :  ndarray of shape (n_samples, 20+lambda_)
        Array containing pseudo amino acid composition.
    
    desc : list of length 20+lambda_
        Order of amino acids and lambda values corresponding to columns in arr.

    References
    ----------

    Tanford C. Contribution of hydrophobic interactions to the stability of the
    globular conformation of proteins. J Am Chem Soc. 1962;84:4240–4247.

    Hopp TP, Woods KR. Prediction of protein antigenic determinants from amino 
    acid sequences. Proc Natl Acad Sci U S A. 1981 Jun;78(6):3824-8. 
    doi: 10.1073/pnas.78.6.3824. PMID: 6167991; PMCID: PMC319665.

    Chou KC. Prediction of protein cellular attributes using pseudo‐amino acid 
    composition. Proteins. 2001;43:246‐255.

    Jain et al. TpPred: A Tool for Hierarchical Prediction of Transport
    Proteins Using Cluster of Neural Networks and Sequence Derived Features. 
    International Journal for Computational Biology (IJCB), 2012;1:28-36.

    Examples
    --------

    >>> from protlearn.features import paac
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> paac_comp, desc = paac(seqs, lambda_=3, remove_zero_cols=True)
    >>> paac_comp
    array([[0.62956037, 0.        , 0.        , 0.62956037, 0.62956037,
            0.        , 0.62956037, 0.62956037, 0.10830885, 0.16327632,
            0.09885446],
           [0.        , 1.49423911, 0.74711955, 0.74711955, 0.74711955,
            0.74711955, 0.74711955, 0.        , 0.04009348, 0.08251201,
            0.13027496]])
    >>> desc
    ['A', 'E', 'G', 'K', 'L', 'P', 'R', 'Y', 'lambda1', 'lambda2', 'lambda3']
    """
    
    # input handling
    X = check_input(X)
    
    # load data
    df = pd.read_csv(PATH+'paac.csv')
    data = np.asarray(df.iloc[:,1:])

    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    desc = [aa for aa in amino_acids]

    for n in range(1, lambda_+1):
        desc.append('lambda' + str(n))
    
    # normalization
    for i in range(data.shape[1]):
        mean = np.mean(data[:,i])
        denom = np.sqrt(sum([(j-mean)**2 for j in data[:,i]])/20)
        data[:,i] = [(j-mean)/denom for j in data[:,i]]

    aa_dict = {amino_acids[i]: data[i,:] for i in range(20)}

    # correlation function
    def corr(aa1, aa2, data):
        return sum([(aa_dict[aa1][i] - aa_dict[aa2][i])**2 \
                    for i in range(data.shape[1])]) / data.shape[1]

    # computing pseudo amino acid composition
    arr = np.zeros((len(X),len(amino_acids)+lambda_))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        theta = []
        for n in range(1, lambda_+1):
            theta.append(
                sum([corr(seq[j], seq[j+n], data) \
                     for j in range(len(seq)-n)])/(len(seq)-n))
        
        cnt = Counter(seq)
        pse_aac = [cnt[aa] / (1+w * sum(theta)) for aa in amino_acids]
        pse_aac = pse_aac + [(w*j) / (1+w * sum(theta)) for j in theta]
        arr[i,:] = pse_aac

    # delete zero columns
    if remove_zero_cols:
        cols_zeros = np.where(~arr.any(axis=0))[0]
        arr = np.delete(arr, cols_zeros, axis=1)
        desc = [i for j, i in enumerate(desc) if j not in cols_zeros]

    return arr, desc