import os
import sys
sys.path.insert(0, os.path.dirname(os.getcwd())+'/protclass')

from preprocessing import txt_to_df


def test_conversion():
    "Test txt_to_df conversion"
    
    # load data
    df = txt_to_df('docs/test_seq.txt', 0)
    
    # test labels and df shape
    assert df.columns[0] == 'Sequence'
    assert df.columns[1] == 'Label'
    assert df.shape == (3, 2)
    
    # test sequences
    assert df['Sequence'][0] == 'AGTYLK'
    assert df['Sequence'][1] == 'VTAYKLLLF'
    assert df['Sequence'][2] == 'ARALTYI'
    
    # test labels
    for i in range(df.shape[0]):
        assert df['Label'][i] == 0