import numpy as np
from features import aaindex1

def test_aaindex1():
    "Test AAIndex1"
    
    # load data
    data = ['AGTYLK', 'VCIMMMPFP', 'LRSAHHN', 'AQEEWD'] 
    
    # get aaindex1
    aaind1, desc = aaindex1(data)
    
    # test shape
    assert aaind1.shape == (4, 553)
    
    # test some indices
    ANDN920101 = np.array([4.3, 4.40555, 4.48714, 4.46]) # index 0
    QIAN880126 = np.array([.01166, -.17111, .05857, -.04333]) # index 277
    KARS160122 = np.array([2.014, 5.48522, 2.789, 1.751]) # index -1
    np.testing.assert_equal(np.round(aaind1[0], 3),\
                            np.round(ANDN920101, 3))
    np.testing.assert_equal(np.round(aaind1[277], 3),\
                            np.round(QIAN880126, 3))
    np.testing.assert_equal(np.round(aaind1[-1], 3),\
                            np.round(KARS160122, 3))
    
    # test standardization (zscore)
    aaind1_z, desc = aaindex1(data, 'zscore')
    # test mean = 0
    for i in range(aaind1_z.shape[0]):
        assert abs(round(aaind1_z[:,1].mean())) == 0
    # test std --> 1
    for i in range(aaind1_z.shape[0]):
        assert round(aaind1_z[:,i].std(), 1) ==\
               round(aaind1_z[:,0].std(), 1)
        
    # test standardization (minmax)
    aaind1_mm, desc = aaindex1(data, 'minmax')
    # test minimum and maximum
    for i in range(aaind1_mm.shape[0]):
        assert round(aaind1_mm[:,i].min()) == 0
        assert round(aaind1_mm[:,i].max()) == 1