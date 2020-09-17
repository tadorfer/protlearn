# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
from Bio import SeqIO

def check_input(X):
    """Check if input has the correct type."""
    # correct fasta extensions for proteins
    ext = ['.fasta', '.faa', '.fa']
    if type(X) == str:
        _, extension = os.path.splitext(X)
        
        # fasta format
        if extension in ext:
            # single fasta sequence
            try:
                X = [str(SeqIO.read(X, 'fasta').seq)]
                
            # multiple fasta sequences
            except:
                X = [str(rec.seq) for rec in SeqIO.parse(X, 'fasta')]
                
        else:
            X = [X] 
    
    elif type(X) == list:
        pass
    
    else:
        raise TypeError('Data must be string or list.')
    
    return X

def check_alpha(X):
    """Check that input is alphabetical."""
    if str.isalpha(X) == False:
        raise ValueError('Data must be alphabetical!')

def check_natural(X):
    """Check that sequence is comprised of only natural amino acids."""
    amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    if set(X).issubset(amino_acids) == False:
        raise ValueError("Data contains sequences with unnatural amino acids. "+ 
                         "Consider running preprocessing.remove_unnatural.")
