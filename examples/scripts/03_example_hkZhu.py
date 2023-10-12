#%% Test the function for HK stacking
import os 
import numpy as np
from rf import read_rf
from rfsed.hkZhu import hk, plothk
import h5py, h5
from obspy import read
current_dir = os.path.dirname(os.path.realpath(__file__))
target_dir2 = os.path.sep.join(current_dir.split(os.path.sep)[:-2])
target_dir1=os.path.sep.join(current_dir.split(os.path.sep)[:-1])
staname='OPLO'
# staname = input("Input Station Name:")
w1, w2, w3 = [0.6, 0.3, 0.1]
savepath = target_dir1 + '/plots/HK_Stacking_Zhu/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
#----------------------------------------------------------
#%% Moho Example 
rffile = target_dir2 + '/data/rfstreams_Moho/rfstreams.h5'
rfst = read_rf(rffile, 'H5')
rfstreams = rfst.select(component='R', station=staname)
K= np.linspace(1.65,1.95,121)
H=np.linspace(20,60,201)
Vp=6.9
Result=hk(rfstreams, layer='Moho', stack=True, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
# Result=hk(rfstreams, layer='Moho', stack=False, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
plothk(HKResult=Result, g = [75.,10., 15., 2.5], rmneg = None, savepath=savepath, format = 'jpg')
#----------------------------------------------------------
#%% Sediment Example
# rffile = target_dir2 + '/data/rfstreams_Sed/rfstreams.h5'
# rfst = read_rf(rffile, 'H5')
# rfstreams = rfst.select(component='R', station=staname)
# K= np.linspace(1.65,2.25,201)
# H=np.linspace(0,10,201)
# Vp=2.5
# Result=hk(rfstreams, layer='Sed', stack=False, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
# # Result=hk(rfstreams, layer='Sed', stack=True, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp)
# plothk(HKResult=Result, g = [75.,10., 15., 2.5], rmneg = None, savepath=savepath, format = 'jpg')
# %%
