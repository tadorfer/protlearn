import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protlearn')
import numpy as np

from preprocessing import txt_to_df
from feature_engineering import aaindex3


def test_aaindex3():
    "Test AAIndex3"
    
    # load data
    df = txt_to_df(path+'/tests/docs/test_seq.txt', 0)
    
    # get aaindex2
    aaind3 = aaindex3(df)
    
    # test shape
    assert aaind3.shape == (4, 43)
    
    # test some triangular indices
    TANS760101 = np.array([-4.72, -5.975, -4.18333, -4.04])
    GODA950101 = np.array([np.nan, -.05, -.1333, -.14])
    ZHAC000106 = np.array([.196, -.34875, .46666, .972])
    np.testing.assert_equal(np.round(aaind3['TANS760101'], 3),\
                            np.round(TANS760101, 3))
    # this column contains NaNs
    assert ('GODA950101' in aaind3) == False
    np.testing.assert_equal(np.round(aaind3['ZHAC000106'], 3),\
                            np.round(ZHAC000106, 3))
    
    # test some square indices
    ZHAC000102 = np.array([-.408, -1.415, .475, 1.532])
    ZHAC000103 = np.array([-.052, -.625, .59166, 1.096])
    ZHAC000105 = np.array([-.242, -.72, .17, .952])
    np.testing.assert_equal(np.round(aaind3['ZHAC000102'], 3),\
                            np.round(ZHAC000102, 3))
    np.testing.assert_equal(np.round(aaind3['ZHAC000103'], 3),\
                            np.round(ZHAC000103, 3))
    np.testing.assert_equal(np.round(aaind3['ZHAC000105'], 3),\
                            np.round(ZHAC000105, 3))
    
    # test standardization (zscore)
    aaind3_z = aaindex3(df, 'zscore')
    # test mean = 0
    for i in range(aaind3_z.shape[0]):
        assert abs(round(aaind3_z.iloc[:,1].mean())) == 0
    # test std --> 1
    for i in range(aaind3_z.shape[0]):
        assert round(aaind3_z.iloc[:,i].std(), 1) ==\
               round(aaind3_z.iloc[:,0].std(), 1)
        
    # test standardization (minmax)
    aaind3_mm = aaindex3(df, 'minmax')
    # test minimum and maximum
    for i in range(aaind3_mm.shape[0]):
        assert round(aaind3_mm.iloc[:,i].min()) == 0
        assert round(aaind3_mm.iloc[:,i].max()) == 1