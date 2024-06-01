import unittest
import warnings
from rfsed import ReverbFilter
from rfsed.ReverbFilter import ResonanceFilt, plotfiltrf
from rfsed.hkYu import hkYu, plothkYu
import numpy as np
from rfsed.util import rfMoho_example, save_plot, save_tests
try:
    import matplotlib
except ImportError:
    Matplotlib = None

class ReverbFilter(unittest.TestCase):
    def test_Resonance_Filt(self):
        rfstream = rfMoho_example()
        preonset=10
        self.assertTrue(ResonanceFilt(rfstream, preonset))

    def test_hkYu(self):
        rfstream = rfMoho_example()
        preonset=10
        FilteredRF= ResonanceFilt(rfstream, preonset)
        hkYu(FltResult=FilteredRF, rayp=0.04, HSubSed=np.linspace(20,60,201), KSubSed=np.linspace(1.65,1.95,121), 
             VpMoho=6.9, VpSed= 2.5, VsSed=1.4, w1SubSed=0.6, w2SubSed=0.3, w3SubSed=0.1)  
        # hkYu(FltResult=FilteredRF, rayp=0.04, HSubSed=np.linspace(20,60,201), KSubSed=np.linspace(1.65,1.95,121), 
        #      VpMoho=6.9, VpSed= 2.5, VsSed=1.4, w1SubSed=0.6, w2SubSed=0.3, w3SubSed=0.1)  
    def test_plothkYu(self):
        rfstream = rfMoho_example()
        preonset=10
        FilteredRF= ResonanceFilt(rfstream, preonset)
        savepath=save_tests()
        HKResults=hkYu(FltResult=FilteredRF, rayp=0.04, HSubSed=np.linspace(20,60,201), KSubSed=np.linspace(1.65,1.95,121), 
                VpMoho=6.9, VpSed= 2.5, VsSed=1.4, w1SubSed=0.6, w2SubSed=0.3, w3SubSed=0.1)  
        plothkYu(hkYuResult=HKResults, savepath=savepath)
        
def suite():
    return unittest.makeSuite(ReverbFilter, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')