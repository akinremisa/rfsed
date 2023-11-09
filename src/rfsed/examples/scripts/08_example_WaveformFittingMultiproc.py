#%% Test the Waveform Fitting Function
import os
import numpy as np
from rfsed.rf import read_rf
from rfsed.WaveformFittingMultiproc import WaveformPara, run_waveformfitting, plotbestmodel
from rfsed.util import rfMoho_example, rfSed_example, save_plot

#----------------------------------------------------------
savedir=save_plot()
savepath = savedir + '/WaveformFittingMultiproc/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)

rfstream= rfMoho_example()
VpSed=3.3
Sedthick=1.8
VpCrust = 6.9
gauparameter=1.25
KMoho= np.linspace(1.65,1.95,121)
HMoho=np.linspace(20,60,201)
# KMoho= np.linspace(1.65,1.95,5)
# HMoho=np.linspace(20,60,5)
rayp=0.04
delta = rfstream[0].stats.delta
wtCorr=0.5
wtRMSE=0.3
wtPG=0.2
format='pdf'
nproc=30
ModelParams=WaveformPara(rfstream, Sedthick, VpSed, VpCrust, rayp, KMoho, HMoho, gauparameter)
Results=run_waveformfitting(nproc, HMoho, ModelParams)
plotbestmodel(Results, ModelParams, wtCorr, wtRMSE, wtPG, savepath, format)