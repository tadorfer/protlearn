import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protlearn')

from ..features import composition


def test_composition():
    "Test sequence compositions"
    
    # load data
    data = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD']
    
    # test relative composition
    comp_rel = composition(data, 'relative')
    
    # test if frequencies add to 1
    for i in range(len(data)):
        assert round(comp_rel.iloc[0,:].sum()) == 1
        
    # test decimal places
    comp_dec2 = composition(data, 'relative', round_fraction=2)
    x = str(comp_dec2['A'][0])
    assert x[::-1].find('.') == 2
    comp_dec2 = composition(data, 'relative', round_fraction=5)
    y = str(comp_dec2['A'][0])
    assert y[::-1].find('.') == 5
    
    # test absolute composition
    comp_abs = composition(data, 'absolute')
    
    # test if frequences == sequence length
    all_lengths = [6, 9, 7, 6]
    for i in range(len(data)):
        assert comp_abs.iloc[i,:].sum() == all_lengths[i]