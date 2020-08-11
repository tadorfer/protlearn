from features import aac

def test_aac():
    "Test sequence compositions"
    
    # load data
    data = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD']
    
    # test relative composition
    comp_rel, aa = aac(data, 'relative')
    
    # test if frequencies add to 1
    for i in range(len(data)):
        assert round(comp_rel[0,:].sum()) == 1
    
    # test absolute composition
    comp_abs, aa = aac(data, 'absolute')
    
    # test if frequences == sequence length
    all_lengths = [6, 9, 7, 6]
    for i in range(len(data)):
        assert comp_abs[i,:].sum() == all_lengths[i]