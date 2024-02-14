import unittest
import warnings
from rfsed.WaveformFittingMultiproc import WaveformPara, run_waveformfitting, plotbestmodel
import numpy as np
from rfsed.util import rfMoho_example, save_tests
try:
    import matplotlib
except ImportError:
    Matplotlib = None

class TestWaveformFittingMultiproc(unittest.TestCase):
    def test_WaveformPara(self):
        rfstream= rfMoho_example()
        VpSed=3.3
        Sedthick=1.8
        VpCrust = 6.9
        gaussian=1.25
        savepath=save_tests()
        rayp=0.04
        KMoho= np.linspace(1.65,1.95,121)
        HMoho=np.linspace(20,60,201)
        wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
        WaveformPara(rfstream, Sedthick, VpSed, VpCrust, rayp, KMoho, HMoho, gaussian)

    def test_run_waveformfitting(self):
            rfstream= rfMoho_example()
            VpSed=3.3
            Sedthick=1.8
            VpCrust = 6.9
            gaussian=1.25
            savepath=save_tests()
            rayp=0.04
            KMoho= np.linspace(1.65,1.95,121)
            HMoho=np.linspace(20,60,201)
            nproc=25
            wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
            ModelParams=WaveformPara(rfstream, Sedthick, VpSed, VpCrust, rayp, KMoho, HMoho, gaussian)
            run_waveformfitting(nproc, HMoho, ModelParams)

    def test_plotbestmodel(self):
            rfstream= rfMoho_example()
            VpSed=3.3
            Sedthick=1.8
            VpCrust = 6.9
            gaussian=1.25
            savepath=save_tests()
            rayp=0.04
            KMoho= np.linspace(1.65,1.95,121)
            HMoho=np.linspace(20,60,201)
            nproc=25
            format='pdf'
            wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
            ModelParams=WaveformPara(rfstream, Sedthick, VpSed, VpCrust, rayp, KMoho, HMoho, gaussian)
            Results=run_waveformfitting(nproc, HMoho, ModelParams)
            plotbestmodel(Results, ModelParams, wtCorr, wtRMSE, wtPG, savepath, format)
        

def suite():
    return unittest.makeSuite(TestWaveformFittingMultiproc, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')