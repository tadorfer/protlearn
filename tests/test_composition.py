import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protclass')

from preprocessing import txt_to_df
from feature_engineering import composition


def test_composition():
    "Test sequence compositions"
    
    # load data
    df = txt_to_df(path+'/tests/docs/test_seq.txt', 0)
    
    # test relative composition
    comp_rel = composition(df, 'relative')
    
    # test if frequencies add to 1
    for i in range(df.shape[0]):
        assert round(comp_rel.iloc[0,:].sum()) == 1
        
    # test decimal places
    comp_dec2 = composition(df, 'relative', 2)
    x = str(comp_dec2['A'][0])
    assert x[::-1].find('.') == 2
    comp_dec2 = composition(df, 'relative', 5)
    y = str(comp_dec2['A'][0])
    assert y[::-1].find('.') == 5
    
    # test absolute composition
    comp_abs = composition(df, 'absolute')
    
    # test if frequences == sequence length
    all_lengths = [6, 9, 7]
    for i in range(df.shape[0]):
        assert comp_abs.iloc[i,:].sum == all_lengths[i]