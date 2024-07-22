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
        preonset=10
        VpSed=3.3
        HSed=1.8
        VpCrust = 6.9
        gaussian=1.25
        rayp=0.04
        # KMoho= np.linspace(1.65,1.95,121)
        # HMoho=np.linspace(20,60,201)
        # Use smaller grid search discretization so the test runs faster
        KMoho= np.linspace(1.65,1.95,20)
        HMoho=np.linspace(25,50,20)
        format='pdf'
        wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
        WaveformFitting(rfstream=rfstream, preonset=preonset, HSed=HSed, VpSed=VpSed, VpCrust=VpCrust, rayp=rayp, KMoho=KMoho, HMoho=HMoho, 
                     gaussian=gaussian, wtCorr=wtCorr, wtRMSE=wtRMSE, wtPG=wtPG, savepath=savepath, format=format)

    def test_PlotWaveformFitting(self):
            savepath=save_tests()
            rfstream= rfMoho_example()
            preonset=10
            VpSed=3.3
            HSed=1.8
            VpCrust = 6.9
            gaussian=1.25
            rayp=0.04
            # KMoho= np.linspace(1.65,1.95,121)
            # HMoho=np.linspace(20,60,201)
            # Use smaller grid search discretization so the test runs faster
            KMoho= np.linspace(1.65,1.95,20)
            HMoho=np.linspace(25,50,20)
            wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
            format='pdf'
            FittingResult=WaveformFitting(rfstream=rfstream, preonset=preonset, HSed=HSed, VpSed=VpSed, VpCrust=VpCrust, rayp=rayp, KMoho=KMoho, HMoho=HMoho, 
                     gaussian=gaussian, wtCorr=wtCorr, wtRMSE=wtRMSE, wtPG=wtPG, savepath=savepath, format=format)
            PlotWaveformFitting(FittingResult, wtCorr, wtRMSE, wtPG, savepath, format)
    

def suite():
    return unittest.makeSuite(TestWaveformFitting, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')