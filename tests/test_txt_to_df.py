import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protclass')

from preprocessing import txt_to_df


def test_conversion():
    "Test txt_to_df conversion"
    
    # load data
    df = txt_to_df(path+'/tests/docs/test_seq.txt', 0)
    
    # test labels and df shape
    assert df.columns[0] == 'Sequence'
    assert df.columns[1] == 'Label'
    assert df.shape == (4, 2)
    
    # test sequences
    assert df['Sequence'][0] == 'AGTYLK'
    assert df['Sequence'][1] == 'VCIMMMPFP'
    assert df['Sequence'][2] == 'LRSAHHN'
    assert df['Sequence'][2] == 'AQEEWD'
    
    # test labels
    for i in range(df.shape[0]):
        assert df['Label'][i] == 0