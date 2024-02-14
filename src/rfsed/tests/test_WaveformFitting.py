import unittest
import warnings
from rfsed.WaveformFitting import WaveformFitting, PlotWaveformFitting 
import numpy as np
from rfsed.util import rfMoho_example, save_tests
try:
    import matplotlib
except ImportError:
    Matplotlib = None

class TestWaveformFitting(unittest.TestCase):
    def test_WaveformFitting(self):
        savepath=save_tests()
        rfstream= rfMoho_example()
        VpSed=2.2
        HSed=0.1
        VpCrust = 6.9
        gaussian=1.25
        rayp=0.04
        KMoho= np.linspace(1.65,1.95,121)
        HMoho=np.linspace(20,60,201)
        format='pdf'
        wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
        WaveformFitting(rfstream, HSed, VpSed, VpCrust, rayp, KMoho, HMoho, 
                     gaussian, wtCorr, wtRMSE, wtPG, savepath, format)

    def test_PlotWaveformFitting(self):
            savepath=save_tests()
            rfstream= rfMoho_example()
            VpSed=2.2
            HSed=0.1
            VpCrust = 6.9
            gaussian=1.25
            rayp=0.04
            KMoho= np.linspace(1.65,1.95,121)
            HMoho=np.linspace(20,60,201)
            wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
            format='pdf'
            FittingResult=WaveformFitting(rfstream, HSed, VpSed, VpCrust, rayp, KMoho, HMoho, 
                        gaussian, wtCorr, wtRMSE, wtPG, savepath, format)
            PlotWaveformFitting(FittingResult, wtCorr, wtRMSE, wtPG, savepath, format)
    

def suite():
    return unittest.makeSuite(TestWaveformFitting, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')