.. _introduction:

Introduction
============

What is protlearn?
------------------

*protlearn* is a feature extraction tool for protein sequences. It is a freely available Python
package that allows the user to efficiently extract amino acid sequence features
from proteins or peptides, which can then be used for a variety of downstream 
machine learning tasks.

Package Contents 
----------------

.. image:: protlearn_summary.png
   :alt: Summary of Modules
   :align: center

*protlearn* is currenty comprised of three main modules:

* preprocessing
* feature extraction
* dimensionality reduction

The :ref:`preprocessing` section allows the user to prepare the datasets according to 
specific needs, such as sequence identity or representation (i.e. conversion 
of alphabetical strings to integers). The :ref:`feature_extraction` section can then be 
used to compute amino acid sequence features from the dataset, such as amino acid  
composition or AAIndex-based physicochemical properties. Finally, 
:ref:`dimensionality_reduction` methods are provided to reduce the dimensionality of the computed 
features, which reduces redundancy and alleviates the computational demand for the classifiers.

Input Formats
-------------

*protlearn* was designed to handle both single and multiple sequence inputs in the
following formats:

* string (single)
* .fasta (single)
* list of strings (multiple)
* .fasta (multiple)

Installing protlearn
--------------------

*protlearn* is currently supported and tested on Python 3.6, 3.7, and 3.8 and can 
be installed and upgraded using the following terminal commands:

.. code::

    pip install protlearn  
    pip install --upragde protlearn 