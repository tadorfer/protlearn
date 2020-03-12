import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protclass')
import numpy as np

from preprocessing import txt_to_df
from feature_engineering import aaindex1


def test_aaindex1():
    "Test AAIndex1"
    
    # load data
    df = txt_to_df(path+'/tests/docs/test_seq.txt', 0)
    
    # get aaindex1
    aaind1 = aaindex1(df)
    
    # test shape
    assert aaind1.shape == (4, 566)
    
    # test some indices
    ANDN920101 = np.array([4.3, 4.40555, 4.48714, 4.46])
    QIAN880126 = np.array([.01166, -.17111, .05857, -.04333])
    KARS160122 = np.array([2.014, 5.48522, 2.789, 1.751])
    np.testing.assert_equal(np.round(aaind1['ANDN920101'], 3),\
                            np.round(ANDN920101, 3))
    np.testing.assert_equal(np.round(aaind1['QIAN880126'], 3),\
                            np.round(QIAN880126, 3))
    np.testing.assert_equal(np.round(aaind1['KARS160122'], 3),\
                            np.round(KARS160122, 3))
    
    # test standardization (zscore)
    aaind1_z = aaindex1(df, 'zscore')
    # test mean = 0
    for i in range(aaind1_z.shape[0]):
        assert abs(round(aaind1_z.iloc[:,1].mean())) == 0
    # test std --> 1
    for i in range(aaind1_z.shape[0]):
        assert round(aaind1_z.iloc[:,i].std(), 1) ==\
               round(aaind1_z.iloc[:,0].std(), 1)
        
    # test standardization (minmax)
    aaind1_mm = aaindex1(df, 'minmax')
    # test minimum and maximum
    for i in range(aaind1_mm.shape[0]):
        assert round(aaind1_mm.iloc[:,i].min()) == 0
        assert round(aaind1_mm.iloc[:,i].max()) == 1