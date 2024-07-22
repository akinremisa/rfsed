import unittest
import warnings
from rfsed import hkZhu
from rfsed.hkZhu import hk, plothk
import numpy as np
from rfsed.util import rfMoho_example, rfSed_example, save_plot, save_tests
try:
    import matplotlib
except ImportError:
    Matplotlib = None

class TesthkZhu(unittest.TestCase):
    def test_hk(self):
        rfstream = rfMoho_example()
        preonset=10
        staname='OPLO'
        rfstreams = rfstream.select(component='R', station=staname)
        K= np.linspace(1.65,1.95,121)
        H=np.linspace(20,60,201)
        Vp=6.9
        w1, w2, w3 =[0.6, 0.3, 0.1]
        self.assertTrue(hkZhu.hk(rfstreams, preonset, layer='Moho', stack=True, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp))
    
    def test_plothk(self):
        rfstream = rfMoho_example()
        staname='OPLO'
        preonset=10
        rfstreams = rfstream.select(component='R', station=staname)
        K= np.linspace(1.65,1.95,121)
        H=np.linspace(20,60,201)
        Vp=6.9
        w1, w2, w3 =[0.6, 0.3, 0.1]
        savepath=save_tests()
        HKResult= hkZhu.hk(rfstreams, preonset, layer='Moho', stack=True, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
        hkZhu.plothk(HKResult=HKResult, g = [75.,10., 15., 2.5], rmneg = None,savepath=savepath, format = 'jpg')

def suite():
    return unittest.makeSuite(TesthkZhu, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')