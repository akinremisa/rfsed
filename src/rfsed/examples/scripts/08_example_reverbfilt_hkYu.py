#%% Test the Resonance filtering and Modified HK Stacking of Yu et al 2015
import os 
import numpy as np
from obspy import read
from rfsed.ReverbFilter import Resonance_Filt, plotfiltrf
from rfsed.util import rfMoho_example, rfSed_example, save_plot

savedir=save_plot()
savepath = savedir + '/Resonance_Filtering/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
rfstream= rfMoho_example()
VpSed=2.1
VpMoho=6.9
SedH= 0.6
VsSed= 0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
gaussalp=1.25
rayp = 0.04
w1Sed, w2Sed, w3Sed = [0.6, 0.3, 0.1]
w1SubSed, w2SubSed, w3SubSed = [0.6, 0.3, 0.1]
#----------------------------------------------------------
FilteredRF= Resonance_Filt(rfstream)
plotfiltrf(FilteredRF, savepath, format = 'jpg')

#%% Test the Modified HK Stacking of Yu et al 2015
# Requires result from Resonance filter 
import os 
import numpy as np
from obspy import read
from rfsed.hkYu import hkYu, plothkYu
from rfsed.util import rfMoho_example, rfSed_example, save_plot

savedir=save_plot()
savepath = savedir + '/HK_Yu_Method/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
rfstream= rfMoho_example()
VpSed=2.1
VpMoho=6.9
SedH= 0.6
VsSed= 0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
gaussalp=1.25
rayp = 0.04
w1Sed, w2Sed, w3Sed = [0.6, 0.3, 0.1]
w1SubSed, w2SubSed, w3SubSed = [0.6, 0.3, 0.1]
#----------------------------------------------------------
FilteredRF= Resonance_Filt(rfstream)
# plotfiltrf(FilteredRF, savepath, format = 'jpg')
HKResults=hkYu(FltResult=FilteredRF, rayp=0.04, HSubSed=np.linspace(20,60,201), KSubSed=np.linspace(1.65,1.95,121), 
               HSed=np.linspace(0,10,201), KSed=np.linspace(1.65,2.25,201), VpMoho=6.9, VpSed= 2.5,  
               w1SubSed=0.6, w2SubSed=0.3, w3SubSed=0.1, w1Sed=0.6, w2Sed=0.3, w3Sed=0.1)
plothkYu(hkYuResult=HKResults, savepath=savepath, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg')