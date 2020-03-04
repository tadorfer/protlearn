import numpy as np
import pandas as pd

def test_encoding_string(padding=True):

    # data
    X = pd.DataFrame(['A', 'C', 'M', 'P'], ['D', 'S', 'V'])

    # list of amino acids for integer encoding
    amino_acids = ['A','C','D','E','F','G','H','I','K','L',
                   'M','N','P','Q','R','S','T','V','W','Y'] 
    int_values = np.arange(1, len(amino_acids)+1)

    enc_list = []
    for seq in X['Sequence']:
        seq_split = [aa for aa in seq]
        seq_trans = [amino_acids.index(i)+1 for i in seq_split]
        enc_list.append(np.asarray(seq_trans))
    enc_arr = np.asarray(enc_list)

    if padding == True:
        # get maximum sequence length
        all_len = [len(i) for i in enc_arr]
        max_len = max(all_len)

        # padding up to max_len with zeros
        enc_arr = [np.pad(enc_arr[i], (0, max_len-len(enc_arr[i])),\
                  'constant') for i in range(enc_arr.shape[0])]
        enc_arr = np.asarray(enc_arr) 

    assert enc_arr.shape[0] == 2
    assert enc_arr.shape[1] == 4