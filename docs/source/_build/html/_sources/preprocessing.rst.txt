.. _preprocessing:

Preprocessing 
=============

Integer-encoding
----------------

The function ``integer_encode`` converts amino acids into their corresponding integers 
based on the specified notation. This is particularly useful for preparing a 
sequence-based model such as a long short-term memory (LSTM) or a gated recurrent 
unit (GRU). *protlearn* provides two different notations for this process:

- 20 standard amino acids (IUPAC standard)
- 26 extended amino acids (IUPAC extended)

Amino acids will be encoded using the integers from 1 onwards, with zeros 
reserved for padding. If ``padding=True``, sequences of unequal lengths will be 
posteriorly padded with zeros to the length of the longest sequence in the 
dataset. 

Examples
########

.. code-block:: python

    from protlearn.preprocessing import integer_encode

    seq = 'ARKLYPGPGEERNK'
    enc, aa = integer_encode(seq)

The output of the above code will be:

.. code-block:: text

    >>> enc
    array([ 1, 15,  9, 10, 20, 13,  6, 13,  6,  4,  4, 15, 12,  9])

    >>> aa
    'ACDEFGHIKLMNPQRSTVWY'

Below is an example using multiple sequences and padding:

.. code-block:: python

    from protlearn.preprocessing import integer_encode

    seqs = ['ARKLY', 'EERNPAA', 'QEPGPGLLLK']
    enc, aa = integer_encode(seqs, padding=True)

The output is as follows:

.. code-block:: text

    >>> enc
    array([[ 1, 15,  9, 10, 20,  0,  0,  0,  0,  0],
           [ 4,  4, 15, 12, 13,  1,  1,  0,  0,  0],
           [14,  4, 13,  6, 13,  6, 10, 10, 10,  9]])

    >>> aa
    'ACDEFGHIKLMNPQRSTVWY'

Remove duplicates
-----------------

The function ``remove_duplicates`` detects and removes duplicate sequences in the dataset.

Example
#######

.. code-block:: python

    from protlearn.preprocessing import remove_duplicates

    seqs = ['ARKLY', 'EERNPAA', 'ARKLY', 'QEPGPGLLLK']
    seqs = remove_duplicates(seqs)

The output of the above code will be:

.. code-block:: text

    >>> seqs
    ['EERNPAA', 'QEPGPGLLLK', 'ARKLY']

Remove unnatural amino acids
----------------------------

The function ``remove_unnatural`` detects and removes sequences containing amino 
acids other than the 20 natural ones.

Example
#######

.. code-block:: python

    from protlearn.preprocessing import remove_unnatural

    seqs = ['ARKLY', 'EERNPJAB', 'QEPGPGLLLK']
    seqs = remove_unnatural(seqs)

The output of the above code will be:

.. code-block:: text

    >>> seqs
    ['ARKLY', 'QEPGPGLLLK']