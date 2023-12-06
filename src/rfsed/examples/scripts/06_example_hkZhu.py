#%% Test the function for HK stacking
import os 
import numpy as np
from rfsed.rf.rfstream import read_rf
from rfsed.hkZhu import hk, plothk
import h5py, h5
from obspy import read
from rfsed.util import rfMoho_example, rfSed_example, save_plot
staname='OPLO'
w1, w2, w3 = [0.6, 0.3, 0.1]
savedir=save_plot()
savepath = savedir + '/HK_Stacking_Zhu/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
#----------------------------------------------------------
#%% Moho Example 
rfst = rfMoho_example()
rfstreams = rfst.select(component='R', station=staname)
K= np.linspace(1.65,1.95,121)
H=np.linspace(20,60,201)
Vp=6.9
Result=hk(rfstreams, layer='Moho', stack=True, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
# Result=hk(rfstreams, layer='Moho', stack=False, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
plothk(HKResult=Result, g = [75.,10., 15., 2.5], rmneg = None, savepath=savepath, format = 'jpg')
#----------------------------------------------------------
#%% Sediment Example
# rfst = rfSed_example()
# rfstreams = rfst.select(component='R', station=staname)
# K= np.linspace(1.65,2.25,201)
# H=np.linspace(0,10,201)
# Vp=2.5
# Result=hk(rfstreams, layer='Sed', stack=False, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
# # Result=hk(rfstreams, layer='Sed', stack=True, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
# plothk(HKResult=Result, g = [75.,10., 15., 2.5], rmneg = None, savepath=savepath, format = 'jpg')
# %%
