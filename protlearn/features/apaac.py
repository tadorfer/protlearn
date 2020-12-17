# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
import pandas as pd
from collections import Counter
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

def apaac(X, *, lambda_=30, w=.05, remove_zero_cols=False, start=1, end=None):
    """Amphiphilic pseudo amino acid composition.

    This feature has the same form as the vanilla amino acid composition, but 
    contains much more information that is related to the sequence order of a 
    protein and the distribution of the hydrophobic and hydrophilic amino acids
    along its chain. For the exact method of calculating the amphiphilic
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

    arr : ndarray of shape (n_samples, 20+2*lambda_)
        Array containing amphiphilic pseudo amino acid composition.
    
    desc : list of length 20+2*lambda_
        Order of amino acids and lambda values corresponding to columns in arr.

    References
    ----------

    Tanford C. Contribution of hydrophobic interactions to the stability of the
    globular conformation of proteins. J Am Chem Soc. 1962;84:4240–4247.

    Hopp TP, Woods KR. Prediction of protein antigenic determinants from amino 
    acid sequences. Proc Natl Acad Sci U S A. 1981 Jun;78(6):3824-8. 
    doi: 10.1073/pnas.78.6.3824. PMID: 6167991; PMCID: PMC319665.

    Chou KC. Using amphiphilic pseudo amino acid composition to predict enzyme 
    subfamily classes. Bioinformatics. 2005 Jan 1;21(1):10-9. 
    doi: 10.1093/bioinformatics/bth466. Epub 2004 Aug 12. PMID: 15308540.

    Jain et al. TpPred: A Tool for Hierarchical Prediction of Transport
    Proteins Using Cluster of Neural Networks and Sequence Derived Features. 
    International Journal for Computational Biology (IJCB), 2012;1:28-36.

    Examples
    --------

    >>> from protlearn.features import apaac
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> apaac_comp, desc = apaac(seqs, lambda_=3, remove_zero_cols=True)
    >>> apaac_comp
    array([[ 1.15705278e+00,  0.00000000e+00,  0.00000000e+00,
             1.15705278e+00,  1.15705278e+00,  0.00000000e+00,
             1.15705278e+00,  1.15705278e+00,  1.38909208e-02,
             3.13346727e-02, -8.12700103e-02, -6.96352986e-02,
            -1.82775274e-05, -5.13547885e-02],
           [ 0.00000000e+00,  1.64484218e+00,  8.22421091e-01,
             8.22421091e-01,  8.22421091e-01,  8.22421091e-01,
             8.22421091e-01,  0.00000000e+00,  4.76440749e-02,
             6.11246926e-02,  1.80650812e-02,  5.07150233e-02,
            -1.93053341e-02,  1.93353709e-02]])
    >>> desc
    ['A',
     'E',
     'G',
     'K',
     'L',
     'P',
     'R',
     'Y',
     'lambda_hphob1',
     'lambda_hphil1',
     'lambda_hphob2',
     'lambda_hphil2',
     'lambda_hphob3',
     'lambda_hphil3']

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
        desc.append('lambda_hphob' + str(n))
        desc.append('lambda_hphil' + str(n))
    
    # normalization
    for i in range(data.shape[1]):
        mean = np.mean(data[:,i])
        denom = np.sqrt(sum([(j-mean)**2 for j in data[:,i]])/20)
        data[:,i] = [(j-mean)/denom for j in data[:,i]]

    aa_dict = {amino_acids[i]: data[i,:] for i in range(20)}

    # computing pseudo amino acid composition
    arr = np.zeros((len(X),len(amino_acids)+2*lambda_))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical  
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information
        tau = []
        for n in range(1, lambda_+1):
            for j in range(2):
                tau.append(
                    sum([aa_dict[seq[k]][j]*aa_dict[seq[k+n]][j] \
                         for k in range(len(seq)-n)])/(len(seq)-n))
        
        cnt = Counter(seq)
        apse_aac = [cnt[aa] / (1+w * sum(tau)) for aa in amino_acids]
        apse_aac = apse_aac + [(w*j) / (1+w * sum(tau)) for j in tau]
        arr[i,:] = apse_aac

    # delete zero columns
    if remove_zero_cols:
        cols_zeros = np.where(~arr.any(axis=0))[0]
        arr = np.delete(arr, cols_zeros, axis=1)
        desc = [i for j, i in enumerate(desc) if j not in cols_zeros]

    return arr, desc