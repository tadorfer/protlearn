# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import os
from Bio import SeqIO

def check_input(X):

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
    "Check that input is alphabetical"
    if str.isalpha(X) == True:
        pass
    else:
        raise ValueError('Data must be alphabetical!')