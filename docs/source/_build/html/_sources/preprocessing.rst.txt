.. _preprocessing:

Preprocessing 
=============

Integer-encoding
----------------

The function ``encode`` converts amino acids into their corresponding integers 
based on the specified notation. This is particularly useful for preparing a 
sequence-based model such as a long short-term memory (LSTM) or a gated recurrent 
unit (GRU). *protlearn* provides three different notations for this process:

- 20 standard amino acids (IUPAC standard)
- 26 extended amino acids (IUPAC extended)
- 20 standard amino acids + X for unknown amino acids

Amino acids will be encoded using the integers from 1 onwards, with zeros 
reserved for padding. If ``padding=True``, sequences of unequal lengths will be 
posteriorly padded with zeros to the length of the longest sequence in the 
dataset. 

Examples
########

.. code-block:: python

    from protlearn.preprocessing import encode

    seq = 'ARKLYPGPGEERNK'
    enc, aa = encode(seq)

The output of the above code will be:

.. code-block:: text

    >>> enc
    array([ 1, 15,  9, 10, 20, 13,  6, 13,  6,  4,  4, 15, 12,  9])

    >>> aa
    'ACDEFGHIKLMNPQRSTVWYBXZJUO'

Below is an example using multiple sequences and padding:

.. code-block:: python

    from protlearn.preprocessing import encode

    seqs = ['ARKLY', 'EERNPAA', 'QEPGPGLLLK']
    enc, aa = encode(seqs, padding=True)

The output is as follows:

.. code-block:: text

    >>> enc
    array([[ 1, 15,  9, 10, 20,  0,  0,  0,  0,  0],
           [ 4,  4, 15, 12, 13,  1,  1,  0,  0,  0],
           [14,  4, 13,  6, 13,  6, 10, 10, 10,  9]])

    >>> aa
    'ACDEFGHIKLMNPQRSTVWYBXZJUO'