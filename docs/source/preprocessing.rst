.. _preprocessing:
.. |br| raw:: html

   <br />

Preprocessing 
=============

Integer-encoding
----------------

.. code-block:: text

    protlearn.preprocessing.integer_encode(X, notation='standard', padding=False)

Encode amino acids as integers.

This function converts amino acids into their corresponding integers 
based on the specified notation, starting at 1. Zeros are reserverd for optional
padding. This is particularly useful for preparing a sequence-based model such 
as a long short-term memory (LSTM) or a gated recurrent unit (GRU). 

Parameters
##########

X: string, fasta, or a list thereof
    Dataset of amino acid sequences.

notation: string, default='standard'
    'standard' : 20 natural amino acids |br|
    'extended' : 20 natural + six additional amino acids: |br| 
    B = aspartic acid or asparagine |br|
    X = unknown or 'other' amino acid |br|
    Z = glutamic acid or glutamine |br|
    J = leucine or isoleucine |br|
    U = selenocysteine |br|
    O = pyrrolysine

padding: bool, default=False
    False : sequences are returned in their original lengths |br|
    True : sequences will be padded with zeros to the length of the longest sequence in the dataset

Returns
#######

enc: ndarray of shape (n_samples,) if padding=False or (n_samples, max_len) if padding=True
    Contains the label-encoded peptide sequences.

amino_acids: amino acid order of enc array
    This serves as a lookup for the encoded sequences.


Examples
########

.. code-block:: python

    >>> from protlearn.preprocessing import integer_encode
    >>> seq = 'ARKLYPGPGEERNK'
    >>> enc, aa = integer_encode(seq)
    >>> enc
    array([ 1, 15,  9, 10, 20, 13,  6, 13,  6,  4,  4, 15, 12,  9])
    >>> aa
    'ACDEFGHIKLMNPQRSTVWY'

Below is an example using multiple sequences and padding. If ``padding=True``, 
sequences of unequal lengths will be posteriorly padded with zeros to the length
of the longest sequence in the dataset. 

.. code-block:: python

    >>> from protlearn.preprocessing import integer_encode
    >>> seqs = ['ARKLY', 'EERNPAA', 'QEPGPGLLLK']
    >>> enc, aa = integer_encode(seqs, padding=True)
    >>> enc
    array([[ 1, 15,  9, 10, 20,  0,  0,  0,  0,  0],
           [ 4,  4, 15, 12, 13,  1,  1,  0,  0,  0],
           [14,  4, 13,  6, 13,  6, 10, 10, 10,  9]])
    >>> aa
    'ACDEFGHIKLMNPQRSTVWY'

Notes
#####

Amino acid sequence used for label-encoding were taken from the official
IUPAC amino acid one-letter notation (Extended IUPAC Protein).

Remove duplicates
-----------------

.. code-block:: text

    protlearn.preprocessing.remove_duplicates(X, verbose=1)

Remove duplicate sequences.

This function detects and removes duplicate sequences from the dataset.

Parameters
##########

X: string, fasta, or a list thereof 
    Dataset of amino acid sequences.

verbose: int, default=1
    0 : no information on duplicates is printed |br|
    1 : prints number of duplicates removed |br|
    2 : prints duplicate sequences and number of times present

Returns
#######

Y: list of length n_samples minus the number of duplicates
    Dataset containing only unique sequences.

Examples
########

.. code-block:: python

    >>> from protlearn.preprocessing import remove_duplicates
    >>> seqs = ['ARKLY', 'EERNPAA', 'ARKLY', 'QEPGPGLLLK']
    >>> seqs = remove_duplicates(seqs)
    >>> seqs
    ['EERNPAA', 'QEPGPGLLLK', 'ARKLY']

Remove unnatural amino acids
----------------------------

.. code-block:: text

    protlearn.preprocessing.remove_unnatural(X)

Remove unnatural sequences.

This function removes sequences containing amino acids other than the 20 natural ones.

Parameters
##########

X: string, fasta, or a list thereof
    Dataset of amino acid sequences.

Returns
########

Y: list of length n_samples minus the number of sequences containing unnatural amino acids
    Dataset containing only sequences comprised of natural amino acids.

Examples
########

.. code-block:: python

    >>> from protlearn.preprocessing import remove_unnatural
    >>> seqs = ['ARKLY', 'EERNPJAB', 'QEPGPGLLLK']
    >>> seqs = remove_unnatural(seqs)
    >>> seqs
    ['ARKLY', 'QEPGPGLLLK']