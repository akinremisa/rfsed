import os 
import scipy
import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rfsed.rf import read_rf, RFStream
from rfsed.rf import get_profile_boxes, iter_event_data, IterMultipleComponents
from rfsed.rf.imaging import plot_profile_map
from rfsed.rf.profile import profile
from tqdm import tqdm
from os.path import exists
from os import mkdir
from rfsed.util import save_calculated_RF, read_raw_waveform_data
#------------------------------------------#
savedir=save_calculated_RF()
rffile = savedir + '00_rf_data.h5'
#%%
eqdata = read_raw_waveform_data()
rfst = RFStream()
for stream3c in tqdm(IterMultipleComponents(eqdata, 'onset', 3)):
    if len(stream3c) != 3:
        continue    
    stream3c.trim2(-10, 40, 'onset')
    stream3c.detrend('linear')
    stream3c.detrend('demean')
    stream3c.taper(type = 'hann', max_percentage=0.07)
    stream3c.filter('bandpass', freqmin = 0.05, 
                  freqmax=1.25, corners=4, zerophase=True)
    stream3c.rf(deconvolve='iterative',  rotate='NE->RT', gauss=np.log(1.25))
    #stream3c.plot()
    rfst.extend(stream3c)
rfst.write(rffile, 'H5')