import os
import sys
path = os.environ.get('TRAVIS_BUILD_DIR')
sys.path.insert(0, path+'/protclass')
import numpy as np

from preprocessing import txt_to_df
from feature_engineering import aaindex2


def test_aaindex2():
    "Test AAIndex2"
    
    # load data
    df = txt_to_df(path+'/tests/docs/test_seq.txt', 0)
    
    # get aaindex2
    aaind2 = aaindex2(df)
    
    # test shape
    assert aaind2.shape == (4, 94)
    
    # test some triangular indices
    ALTS910101 = np.array([-2, -.125, .333, -2])
    VOGG950101 = np.array([4.28, 5.2, 6.05, 4.32])
    CROG050101 = np.array([-1.8, .625, .5, -.4])
    np.testing.assert_equal(np.round(aaind2['ALTS910101'], 3),\
                            np.round(ALTS910101, 3))
    np.testing.assert_equal(np.round(aaind2['VOGG950101'], 3),\
                            np.round(VOGG950101, 3))
    np.testing.assert_equal(np.round(aaind2['CROG050101'], 3),\
                            np.round(CROG050101, 3))
    
    # test some square indices
    LINK010101 = np.array([.0266, .0955, .13, .1276])
    KOSJ950108 = np.array([1.62, 18.0875, 16.05, 14.68])
    DOSZ010101 = np.array([1.32, 15.7625, -1.1833, -5.12])
    np.testing.assert_equal(np.round(aaind2['LINK010101'], 3),\
                            np.round(LINK010101, 3))
    np.testing.assert_equal(np.round(aaind2['KOSJ950108'], 3),\
                            np.round(KOSJ950108, 3))
    np.testing.assert_equal(np.round(aaind2['DOSZ010101'], 3),\
                            np.round(DOSZ010101, 3))
    
    # test standardization (zscore)
    aaind2_z = aaindex2(df, 'zscore')
    # test mean = 0
    for i in range(aaind2_z.shape[0]):
        assert abs(round(aaind2_z.iloc[:,1].mean())) == 0
    # test std --> 1
    for i in range(aaind2_z.shape[0]):
        assert round(aaind2_z.iloc[:,i].std(), 1) ==\
               round(aaind2_z.iloc[:,0].std(), 1)
        
    # test standardization (minmax)
    aaind2_mm = aaindex2(df, 'minmax')
    # test minimum and maximum
    for i in range(aaind2_mm.shape[0]):
        assert round(aaind2_mm.iloc[:,i].min()) == 0
        assert round(aaind2_mm.iloc[:,i].max()) == 1