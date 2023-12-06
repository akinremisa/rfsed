#%% Test the Waveform Fitting Function
import os 
from os.path import exists
from os import mkdir
import numpy as np
from glob import glob
from rfsed.rf import read_rf, RFStream
from rfsed.WaveformFitting import WaveformFitting, PlotWaveformFitting    
from rfsed.util import rfMoho_example, rfSed_example, save_plot

#----------------------------------------------------------
savedir=save_plot()
savepath = savedir + '/WaveformFitting/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)

rfstream= rfMoho_example()
VpSed=2.2
HSed=0.1
VpCrust = 6.9
gaussian=1.25
rayp=0.04
KMoho= np.linspace(1.65,1.95,121)
HMoho=np.linspace(20,60,201)
# KMoho= np.linspace(1.65,1.95,5)
# HMoho=np.linspace(20,60,5)
wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
format='pdf'
FittingResult = WaveformFitting(rfstream, HSed, VpSed, VpCrust, rayp, KMoho, HMoho, 
                     gaussian, wtCorr, wtRMSE, wtPG, savepath, format)
#%%
PlotWaveformFitting(FittingResult, wtCorr, wtRMSE, wtPG, savepath, format)
#%%