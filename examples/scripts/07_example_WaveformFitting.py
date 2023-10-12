#%% Test the Waveform Fitting Function
import os 
from os.path import exists
from os import mkdir
import numpy as np
from glob import glob
from rf import read_rf, RFStream
from rfsed.WaveformFitting import WaveformFitting, PlotWaveformFitting    
current_dir = os.path.dirname(os.path.realpath(__file__))
target_dir2 = os.path.sep.join(current_dir.split(os.path.sep)[:-2])
target_dir1=os.path.sep.join(current_dir.split(os.path.sep)[:-1])
#----------------------------------------------------------
savepath = target_dir1 + '/plots/WaveformFitting/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
rffileMoho = target_dir2 + '/data/rfstreams_Moho/rfstreams.h5'
rfstream= read_rf(rffileMoho, 'H5')
VpSed=2.2
HSed=0.1
VpCrust = 6.9
gaussian=1.25
rayp=0.04
# KMoho= np.linspace(1.65,1.95,121)
# HMoho=np.linspace(20,60,201)
KMoho= np.linspace(1.65,1.95,5)
HMoho=np.linspace(20,60,5)
wtCorr, wtRMSE, wtPG = [0.5, 0.3, 0.2] #wtCorr+wtRMSE+wtPG=1.0
format='pdf'
FittingResult = WaveformFitting(rfstream, HSed, VpSed, VpCrust, rayp, KMoho, HMoho, 
                     gaussian, wtCorr, wtRMSE, wtPG, savepath, format)
#%%
PlotWaveformFitting(FittingResult, wtCorr, wtRMSE, wtPG, savepath, format)
#%%