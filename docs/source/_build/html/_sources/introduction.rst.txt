.. _introduction:

Introduction
============

What is ProtLearn?
------------------

ProtLearn is a protein feature extraction tool for machine learning tasks. It 
is a freely available Python package that allows the user to efficiently extract
amino acid features from proteins or peptides, which can then be used for a 
variety of downstream tasks such as the prediction of protein families or 
secondary and tertiary structures.

Package Contents 
----------------

ProtLearn's main contents are divided into three main sections:

* Preprocessing
* Feature extraction
* Dimensionality Reduction

The :ref:`preprocessing` section allows the user to prepare the datasets according to 
specific needs, such as sequence redundancy or representation (i.e. conversion 
of alphabetical strings to integers). The :ref:`feature_extraction` section can then be 
used to compute amino acid features from the dataset, such as amino acid  
composition or AAIndex-based physicochemical properties. Finally, 
:ref:`dimensionality_reduction` methods are provided to reduce the dimensionality of the computed 
features, which alleviates the computational demand for the classifiers.

Input Formats
-------------

ProtLearn was designed to handle both single and multiple sequence inputs in the
following formats:

* string (single)
* .fasta (single)
* list of strings (multiple)
* .fasta (multiple)

Installing ProtLearn
--------------------

ProtLearn is currently supported and tested on Python 3.6, 3.7, and 3.8 and can 
be installed and upgraded using the following terminal commands:

.. code::

    pip install protlearn  
    pip install --upragde protlearn 