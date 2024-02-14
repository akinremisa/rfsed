import unittest
import warnings
from rfsed import hkSeqYeck
from rfsed.hkSeqYeck import hkSeq, plotSeqhk
import numpy as np
from rfsed.util import rfMoho_example, rfSed_example, save_plot, save_tests
try:
    import matplotlib
except ImportError:
    Matplotlib = None

class hkSeqYeck(unittest.TestCase):
    def test_hkSeq(self):
        rfstSed = rfSed_example()
        preonset=10
        staname='OPLO'
        rfstreamSed = rfstSed.select(component='R', station=staname)
        KSed= np.linspace(1.65,2.25,201)
        HSed=np.linspace(0,10,201)
        VpSed=2.5
        w1Sed, w2Sed, w3Sed = [0.6, 0.3, 0.1]
        rfstMoho = rfMoho_example()
        rfstreamMoho = rfstMoho.select(component='R', station=staname)
        KMoho= np.linspace(1.65,1.95,121)
        HMoho=np.linspace(20,60,201)
        VpMoho=6.9
        w1Moho, w2Moho, w3Moho = [0.6, 0.3, 0.1]
        self.assertTrue(hkSeq(rfstreamSed, rfstreamMoho,  preonset, w1Sed = w1Sed, 
             w2Sed = w2Sed, w3Sed=w3Sed, KSed=KSed, HSed=HSed, VpSed=VpSed,
             w1Moho = w1Moho, w2Moho = w2Moho, w3Moho=w3Moho, KMoho=KMoho, 
             HMoho=HMoho, VpMoho=VpMoho, stack = False))
    def test_plotSeqhk(self):
        rfstSed = rfSed_example()
        preonset=10
        staname='OPLO'
        rfstreamSed = rfstSed.select(component='R', station=staname)
        KSed= np.linspace(1.65,2.25,201)
        HSed=np.linspace(0,10,201)
        VpSed=2.5
        w1Sed, w2Sed, w3Sed = [0.6, 0.3, 0.1]
        rfstMoho = rfMoho_example()
        rfstreamMoho = rfstMoho.select(component='R', station=staname)
        KMoho= np.linspace(1.65,1.95,121)
        HMoho=np.linspace(20,60,201)
        VpMoho=6.9
        w1Moho, w2Moho, w3Moho = [0.6, 0.3, 0.1]
        savepath=save_tests()
        SequentialHKResult=hkSeq(rfstreamSed, rfstreamMoho, preonset, w1Sed = w1Sed, 
             w2Sed = w2Sed, w3Sed=w3Sed, KSed=KSed, HSed=HSed, VpSed=VpSed,
             w1Moho = w1Moho, w2Moho = w2Moho, w3Moho=w3Moho, KMoho=KMoho, 
             HMoho=HMoho, VpMoho=VpMoho, stack = False)
        plotSeqhk(SequentialHKResult=SequentialHKResult, savepath=savepath)
        

def suite():
    return unittest.makeSuite(hkSeqYeck, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')