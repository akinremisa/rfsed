import unittest
import warnings
from rfsed import ReverbFilter
from rfsed.ReverbFilter import Resonance_Filt, plotfiltrf
import numpy as np
from rfsed.util import rfMoho_example, save_plot, save_tests
try:
    import matplotlib
except ImportError:
    Matplotlib = None

class ReverbFilter(unittest.TestCase):
    def test_Resonance_Filt(self):
        rfstream = rfMoho_example()
        self.assertTrue(Resonance_Filt(rfstream))
    def test_plotfiltrf(self):
        rfstream = rfMoho_example()
        savepath=save_tests()
        FilteredRF= Resonance_Filt(rfstream)
        plotfiltrf(FilteredRF, savepath, format = 'jpg')        

def suite():
    return unittest.makeSuite(ReverbFilter, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')