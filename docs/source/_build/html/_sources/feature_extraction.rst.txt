.. _feature_extraction:
.. |br| raw:: html

   <br />

Feature Extraction 
==================

aac 
---

.. code-block:: text

    protlearn.features.aac(X, method='relative', start=1, end=None)

Amino acid composition.

This function returns the frequency of amino acids for each sequence in the dataset. 

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

aac :  ndarray of shape (n_samples, n_unique_amino_acids)
    Array containing the amino acid composition.

amino_acids : amino acid order of aac array
    Corresponds to the columns in aac.

Examples
########

.. code-block:: python

    >>> from protlearn.features import aac
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> comp, aa = integer_encode(seqs)
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

Coming soon!

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

entropy
-------

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

motif 
-----

Coming soon!

ngram 
-----

Coming soon!

paac 
----

Coming soon!

posrich 
-------

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