.. _feature_extraction:
.. |br| raw:: html

   <br />

Feature Extraction 
==================

length 
------

.. code-block:: text

    protlearn.features.length(X, *, method='int')

Sequence length in amino acids.

The number of amino acids that a protein or peptide is comprised of will be 
counted and returned as either an integer (single sequence), array of 
integers, or a one-hot-encoded array.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.

method: string, default='int'
    'int' : interval data |br|
    'ohe' : one-hot encoded data

Returns
#######

arr: int or ndarray of shape (n_samples, ) or (n_samples, n_unique_lengths) 
    Array containing sequence lengths.

Examples
########

.. code-block:: python

    >>> from protlearn.features import length
    >>> seqs = ['ARKLY', 'EERKPGL', 'LLYPGP']
    >>> l_int = length(seqs)
    >>> l_int
    array([5, 7, 6])
    >>> l_ohe = length(seqs, method='ohe')
    array([[1., 0., 0.],
           [0., 0., 1.],
           [0., 1., 0.]])

aac 
---

.. code-block:: text

    protlearn.features.aac(X, *, method='relative', start=1, end=None)

Amino acid composition.

This function returns the frequency of amino acids for each sequence in the dataset. 

`aac` is calculated as follows:

.. math::

   aac(i) = \frac{ n_i }{N}

where *i* denotes the 20 amino acid residues, *n*\ :sub:`i` \ is the frequency of each 
residue i, and *N* is the total number of residues in the sequence.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.

method: string, default='relative'
    'absolute' : absolute amino acid composition |br|
    'relative' : relative amino acid composition

start: int, default=1
    Determines the starting point of the amino acid sequence. This number is based on one-based indexing.

end: int, default=None
    Determines the end point of the amino acid sequence. Similarly to start, this number is based on one-based indexing.


Returns
#######

arr:  ndarray of shape (n_samples, n_unique_amino_acids)
    Array containing the amino acid composition.

amino_acids: string
    Corresponds to the columns (amino acids) in arr.

Examples
########

.. code-block:: python

    >>> from protlearn.features import aac
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> comp, aa = aac(seqs)
    >>> comp
    array([[0.2       , 0.        , 0.        , 0.2       , 0.2       ,
            0.        , 0.2       , 0.2       ],
           [0.        , 0.28571429, 0.14285714, 0.14285714, 0.14285714,
            0.14285714, 0.14285714, 0.        ]])
    >>> aa
    'AEGKLPRY'

Note that columns containing all zeros have been removed from the final array.

aaindex1
--------

.. code-block:: text

    protlearn.features.aaindex1(X, *, standardize='none', start=1, end=None)

AAIndex1-based physicochemical properties.

AAindex1 ver.9.2 (release Feb, 2017) is a set of 20 numerical values 
representing various physicochemical and biological properties of amino 
acids. Currently, it contains 566 indices, of which 553 contain no NaNs. 
The indices will be collected for each amino acid in the sequence, 
then averaged across the sequence. 

`aaindex1` is calculated as follows:

.. math::

   aaindex1(i) = \frac{ \sum_{n=1}^{N}AAindex_i (aa_n) }{N}

where *i* denotes the 566 AAIndex1 indices, *aa*\ :sub:`n` \ denotes the amino acid at 
position *n*, and N is the total number of residues in the sequence.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.

standardize: string, default='none'
    'none' : unstandardized index matrix will be returned |br|
    'zscore' : index matrix is standardized to have a mean of 0 and standard deviation of 1. |br|
    'minmax' : index matrix is normalized to have a range of [0, 1].

start: int, default=1
    Determines the starting point of the amino acid sequence. This number is based on one-based indexing.

end: int, default=None
    Determines the end point of the amino acid sequence. Similarly to start, this number is based on one-based indexing.

Returns
#######

arr: ndarray of shape (n_samples, 553-566) 
    Array containing the AAIndex1 physicochemical properties.

desc: list of length 553-566
    Corresponds to the columns (AAIndices) in arr.

Examples
########

.. code-block:: python

    >>> from protlearn.features import aaindex1
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> aaind, inds = aaindex1(seqs, standardize='zscore')
    >>> aaind.shape
    (2, 553)
    >>> len(inds)
    553

Notes
#####

Columns (indices) containing NaNs will be removed. Thus, the resulting index
matrix will have a column size between 553-566.

References
##########

- Nakai, K., Kidera, A., and Kanehisa, M.; Cluster analysis of amino acid indices for prediction of protein structure and function. Protein Eng. 2, 93-100 (1988). [PMID:3244698]
- Tomii, K. and Kanehisa, M.; Analysis of amino acid indices and mutation matrices for sequence comparison and structure prediction of proteins. Protein Eng. 9, 27-36 (1996). [PMID:9053899]
- Kawashima, S., Ogata, H., and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 27, 368-369 (1999). [PMID:9847231]
- Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 28, 374 (2000). [PMID:10592278]
- Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., and Kanehisa, M.; AAindex: amino acid index database, progress report 2008. Nucleic Acids Res. 36, D202-D205 (2008). [PMID:17998252]

ngram
-----

.. code-block:: text

    protlearn.features.ngram(X, *, method='relative', start=1, end=None)

N-gram composition.

This function computes the di- or tripeptide composition of amino acid 
sequences. Therefore, the function parameter *n* can only take on 
the arguments 2 and 3 - otherwise, it will raise a ValueError.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.
    
n: int, default=2
    Integer denoting the desired n-gram composition. |br|
    2 : dipeptide composition |br|
    3 : tripepitde composition
    
method: string, default='relative'
    'absolute': absolute n-gram composition |br|
    'relative': relative n-gram composition

start: int, default=1
    Determines the starting point of the amino acid sequence. This number is
    based on one-based indexing.

end: int, default=None
    Determines the end point of the amino acid sequence. Similarly to start,
    this number is based on one-based indexing.
    
Returns
#######

arr: ndarray of shape (n_samples, n_unique^n)
    Depending on n, the returned array will be of size: |br|
    - (n_samples, 400) for dipeptide composition |br|
    - (n_samples, 8000) for tripeptide composition |br|
    if all possible n-gram combinations are represented.

n-grams: list of length n_unique^n
    List of n-grams corresponding to columns in arr.

Examples
########

.. code-block:: python

    >>> from protlearn.features import ngram
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> di, ngrams = ngram(seqs, n=2)
    >>> di
    array([[0.25      , 0.25      , 0.25      , 0.25      , 0.        ,
            0.        , 0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.16666667, 0.16666667,
            0.16666667, 0.16666667, 0.16666667, 0.16666667]])
    >>> ngrams
    ['AR', 'KL', 'LY', 'RK', 'EE', 'ER', 'GL', 'KP', 'PG']
    >>> tri, ngrams = ngram(seqs, n=3)
    array([[0.33333333, 0.33333333, 0.33333333, 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.2       , 0.2       ,
            0.2       , 0.2       , 0.2       ]])
    >>> ngrams
    ['ARK', 'KLY', 'RKL', 'EER', 'ERK', 'KPG', 'PGL', 'RKP']

entropy
-------

.. code-block:: text

    protlearn.features.entropy(X, *, standardize='none', start=1, end=None)

Shannon entropy.

This function computes the Shannon entropy for each sequence in the 
dataset as follows:

.. math::

   H(X) = -\sum_{i=1}^{20}P(x_i)log_2 P(x_i)

where *i* denotes the 20 amino acids and *P(x*\ :sub:`i`\) denotes the 
probability of a given amino acid in the sequence.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.

standardize: string, default='none'
    'none' : unstandardized index matrix will be returned |br|
    'zscore' : index matrix is standardized to have a mean of 0 and standard deviation of 1. |br|
    'minmax' : index matrix is normalized to have a range of [0, 1].

start: int, default=1
    Determines the starting point of the amino acid sequence. This number is
    based on one-based indexing.

end: int, default=None
    Determines the end point of the amino acid sequence. Similarly to start,
    this number is based on one-based indexing.

Returns
#######

arr:  ndarray of shape (n_samples,) if len(X) > 1, otherwise float
    Array containing Shannon entropy values for each sequence.

Examples
########

.. code-block:: python

    >>> from protlearn.features import entropy
    >>> seqs = ['ARKLY', 'EERKPGL', 'AAAAAALY']
    >>> ent = entropy(seqs)
    >>> ent
    array([2.32192809, 2.52164064, 0.64020643])

posrich
-------

.. code-block:: text

    protlearn.features.posrich(X, *, position, aminoacid)

Position-specific amino acids.

This function returns a binary vector or matrix in which ones indicate the 
presence of the given amino acid(s) at the specified position(s), and zeros 
indicate their absence. 

Parameters
##########

X: string, fasta, or a list thereof
    Dataset of amino acid sequences.
    
position: int or list
    Integer or list of integers denoting the position(s) in the sequence. 

aminoacid: string or list
    String or list of strings indicating the amino acid(s) of interest.
    
Returns
#######

arr: ndarray of shape (n_samples, ) or (n_samples, n_positions)
    Binary vector indicating position-specific presence of amino acids.

Examples
########

.. code-block:: python

    >>> from protlearn.features import posrich
    >>> seqs = ['ARKLY', 'ERNLAPG', 'YRLQLLLY']   
    >>> pos_single = posrich(seqs, position=4, aminoacid='L')
    >>> pos_single
    array([1., 1., 0.])
    >>> pos_multiple = posrich(seqs, position=[2,3,4], aminoacid=['R','N','L'])
    array([[1., 0., 1.],
           [1., 1., 1.],
           [1., 0., 0.]])

Notes
#####

The position argument is based on one-based indexing.

motif
-----

.. code-block:: text

    protlearn.features.motif(X, pattern, *, start=1, end=None)

Sequence motifs.

This function returns a binary vector indicating the presence of a specified 
amino acid sequence motif.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.
    
pattern: string
    Represents the sequence motif. |br|
    x --> any amino acid |br|
    [XY] --> X or Y |br|
    {X} --> any amino acid except X

start: int, default=1
    Determines the starting point of the amino acid sequence. This number is
    based on one-based indexing.

end: int, default=None
    Determines the end point of the amino acid sequence. Similarly to start,
    this number is based on one-based indexing.

Returns
#######

arr:  ndarray of shape (n_samples,)
    Binary vector indicating the presence of the motif in sequences.

Examples
########

.. code-block:: python

    >>> from protlearn.features import motif
    >>> seqs = ['AARKYLL', 'LELCDPGPG', 'RAAANCDD']  
    >>> pattern1 = pattern = 'AAx[KC]'
    >>> m1 = motif(seqs, pattern1)
    >>> m1
    array([1., 0., 1.])
    >>> pattern2 = 'xxC[DA]xx{Y}'
    >>> m2 = motif(seqs, pattern2)
    >>> m2
    array([0., 1., 0.])    

Notes
#####

Based on the example above, 'pattern1' is interpreted as follows:
Two consecutive amino acids 'A', followed by any amino acid, followed by
either a 'K' or a 'C'. 

Likewise, pattern2 is interpreted as follows:
Any two consecutive amino acids, followed by a 'C', followed by either a 'D'
or an 'A', followed by any two amino acids, followed by any amino acid
except 'Y'.

apaac
-----

Coming soon!

atc 
---

Coming soon!

binary 
------

Coming soon!

csksaap 
-------

Coming soon!

ctd 
---

Coming soon!

ctdc
----

Coming soon!

ctdd
----

Coming soon!

ctdt
----

Coming soon!

geary 
-----

Coming soon!

moran 
-----

Coming soon!

moreau_broto
------------

Coming soon!

paac 
----

Coming soon!

qso 
---

Coming soon!

socn 
----

Coming soon!

.. 
    Pseudo amino acid composition
    # Chou, K.C., 2001, Prediction of Protein Cellular Attributes Using PseudoAmino Acid Composition, Proteins: Structure, Function, and Genetics
    # Hydrophobicity values: Tanford C., Contribution of Hydrophobic Interactions to the Stability of the Globular Conformation of Proteins, J. Am. Chem. Soc. 84:4240-4274(1962)
    # Hydrophilicity values: Hopp & Woods, Prediction of protein antigenic determinants from amino acid sequences, PNAS 1981
    # Side chain mass: Jain et al., 2012, TpPred: A Tool for Hierarchical Prediction of Transport Proteins Using Cluster of Neural Networks and Sequence Derived Features

.. Conjoint triad descriptors
    # Predicting protein-protein interactions based only on sequences information
    # PyBioMed: a python library for various molecular representations of chemicals, proteins and DNAs and their interactions

.. Atomic and bond composition
    # An in silico platform for predicting, screening and designing of antihypertensive peptides

.. Binary profile pattern
    # Identification of conformational B-cell Epitopes in an antigen from its primary sequence

.. Moran Autocorrelation
    # protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences (default features 8)