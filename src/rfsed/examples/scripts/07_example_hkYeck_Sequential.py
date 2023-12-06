#%% Test the function Sequential Stacking for H-K
import os 
import numpy as np
from rfsed.rf.rfstream import read_rf
from rfsed.hkSeqYeck import hkSeq, plotSeqhk
from rfsed.util import rfMoho_example, rfSed_example, save_plot

savedir=save_plot()
savepath = savedir + '/HK_Sequential_Yeck/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
staname='OPLO'
#----------------------------------------------------------
# Sediment Parameters
rfstSed = rfSed_example()
rfstreamSed = rfstSed.select(component='R', station=staname)
KSed= np.linspace(1.65,2.25,201)
HSed=np.linspace(0,10,201)
VpSed=2.5
w1Sed, w2Sed, w3Sed = [0.6, 0.3, 0.1]
#----------------------------------------------------------
# Moho Parameters
rfstMoho = rfMoho_example()
rfstreamMoho = rfstMoho.select(component='R', station=staname)
KMoho= np.linspace(1.65,1.95,121)
HMoho=np.linspace(20,60,201)
VpMoho=6.9
w1Moho, w2Moho, w3Moho = [0.6, 0.3, 0.1]
#----------------------------------------------------------
SequentialHKResult = hkSeq(rfstreamSed, rfstreamMoho,  w1Sed = w1Sed, 
             w2Sed = w2Sed, w3Sed=w3Sed, KSed=KSed, HSed=HSed, VpSed=VpSed,
             w1Moho = w1Moho, w2Moho = w2Moho, w3Moho=w3Moho, KMoho=KMoho, 
             HMoho=HMoho, VpMoho=VpMoho, stack = False)
plotSeqhk(SequentialHKResult=SequentialHKResult, g = [75.,10., 15., 2.5], 
          rmneg = None, savepath=savepath, format = 'jpg')
#%% Stacked RF Example
# SequentialHKResult = hkSeq(rfstreamSed, rfstreamMoho,  w1Sed = w1Sed, 
#              w2Sed = w2Sed, w3Sed=w3Sed, KSed=KSed, HSed=HSed, VpSed=VpSed,
#              w1Moho = w1Moho, w2Moho = w2Moho, w3Moho=w3Moho, KMoho=KMoho, 
#              HMoho=HMoho, VpMoho=VpMoho, stack = True)
# plotSeqhk(SequentialHKResult=SequentialHKResult, g = [75.,10., 15., 2.5], 
#           rmneg = None, savepath=savepath, format = 'jpg')