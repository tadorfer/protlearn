# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
from collections import Counter
from ..utils.validation import check_input

def ctd(X, start=1, end=None):
    """Compute conjoint triad descriptors.

        Parameters
        ----------

        X : string, fasta, or a list thereof 

        start : int, default=1
            Determines the starting point of the amino acid sequence.

        end : int, default=None
            Determines the end point of the amino acid sequence.

        Returns
        -------

        arr :  ndarray of shape (n_samples, 343)

        ctd_list : list of triad classes

        """

    # input handling
    X = check_input(X)

    # define classes
    classes = {'A': 1, 'G': 1, 'V': 1,
               'I': 2, 'L': 2, 'F': 2, 'P': 2,
               'Y': 3, 'M': 3, 'T': 3, 'S': 3,
               'H': 4, 'N': 4, 'Q': 4, 'W': 4,
               'R': 5, 'K': 5,
               'D': 6, 'E': 6,
               'C': 7}

    # compute CTD
    ctd = dict()
    for i, seq in enumerate(X):
        if str.isalpha(seq) == True:
            pass
        else:
            raise ValueError('Data must be alphabetical!')
        seq = seq[start-1:end]
        seq = ''.join([str(classes[aa]) for aa in seq])
        keys = [seq[x:x+3] for x in range(len(seq)-2)]
        unq = sorted(set(keys))
        ctd_vals = sorted(Counter(keys).items())
        vals = [i[1] for i in ctd_vals]

        for num, j in enumerate(unq):
            if j in ctd:
                ctd[j].append(vals[num])
            else:
                ctd[j] = i*[0]+[vals[num]]

        # append values not present in ctd with zero
        if i != 0:
            maxlen = max([len(c) for c in ctd.values()])
            for z in ctd.values():
                if len(z) < maxlen:
                    z.append(0)
    
    arr = np.array(list(ctd.values()), dtype=float).T
    ctd_list = list(ctd.keys())
    
    return arr, ctd_list