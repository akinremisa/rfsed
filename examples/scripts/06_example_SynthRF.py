#%% Importing the libraries
import numpy as np
import rffw 
import os 
from rfsed.SynthRF import hkSynth, plothkSynth, hkSeqSynth, plothkSeqSynth, ResonanceFilt, plotfiltSynthrf
current_dir = os.path.dirname(os.path.realpath(__file__))
target_dir2 = os.path.sep.join(current_dir.split(os.path.sep)[:-2])
target_dir1=os.path.sep.join(current_dir.split(os.path.sep)[:-1])
#%% Test the Synthetic HK function
#----------------------------------------------------------
# Synthetic Parameters
moddim = 4
rayp = 0.05826535137278129   #for 60 degree distance
gaussalp = 1.25
delay = 10
n = 2001
delta = 0.025
#----------------------------------------------------------
# Model Original
#----------------------------------------------------------
depth = np.array([35, 77.5])
vp = np.array([6.90, 8.045])
vs = np.array([3.951172240000002, 4.49])
rho = np.array([2.939564099940003, 3.299])
rf= rffw.rffw(depth, vp, rho, vs, rayp, gaussalp, delay, n, delta)
l = len(rf)
t = np.arange(0, l)
t = (delta *  t) - delay
Vp=vp[0]
K= np.linspace(1.65,1.95,121)
H=np.linspace(20,60,201)
w1, w2, w3 = [0.6, 0.3, 0.1]
savepath = target_dir1 + '/plots/SynthRF_HK_Result/'
if not os.path.exists(savepath):  # create data folder if necessary
    os.mkdir(savepath)
#-------------------------------------------------------------------------
HKResultSynth=hkSynth(Synthrf=rf, time=t, rayp=rayp, w1 = w1, w2 = w2, w3 = w3, K= K, H=H, Vp = Vp, layer = None)
plothkSynth(HKResultSynth, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg', savepath=savepath)

#%% Test the Synthetic Sequential HK function
# import numpy as np
# import rffw
# import os 
# #----------------------------------------------------------
# # Synthetic Parameters
# moddim = 4
# rayp = 0.05826535137278129   #for 60 degree distance
# gaussalp_Moho = 1.25
# gaussalp_Sed = 5
# delay = 10
# n = 2001
# delta = 0.025
# #----------------------------------------------------------
# # Model Original
# #----------------------------------------------------------
# depth = np.array([2, 35, 77.5])
# vp = np.array([3, 6.90, 8.045])
# VpSed=vp[0]
# VpCrust=vp[1]
# VpMoho=vp[1]
# SedH=depth[0]
# VsSed=0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
# SedDen=((1.6612*(VpSed)) - (0.4721 * ((VpSed)**2)) + (0.0671* ((VpSed)**3)) - (0.0043* ((VpSed)**4)) + (0.000106* ((VpSed)**5)))
# VsCrust=0.7858 - 1.2344*VpCrust + 0.7949*VpCrust**2 - 0.1238*VpCrust**3 + 0.0064*VpCrust**4
# CrustDen=((1.6612*(VpCrust)) - (0.4721 * ((VpCrust)**2)) + (0.0671* ((VpCrust)**3)) - (0.0043* ((VpCrust)**4)) + (0.000106* ((VpCrust)**5)))
# vs = np.array([VsSed,VsCrust, 4.49])
# rho = np.array([SedDen, CrustDen, 3.299])
# rfSynthMoho= rffw.rffw(depth, vp, rho, vs, rayp, gaussalp_Moho, delay, n, delta)
# rfSynthSed= rffw.rffw(depth, vp, rho, vs, rayp, gaussalp_Sed, delay, n, delta)
# l = len(rfSynthMoho)
# t = np.arange(0, l)
# t = (delta *  t) - delay
# w1, w2, w3 = [0.6, 0.3, 0.1]
# savepath = target_dir1 + '/plots/SynthRF_SeqHK_Result/'
# if not os.path.exists(savepath):  # create data folder if necessary
#     os.mkdir(savepath)
# KSed= np.linspace(1.65,2.25,201)
# HSed=np.linspace(0,10,201)
# w1Sed, w2Sed, w3Sed = [0.6, 0.3, 0.1]
# KMoho= np.linspace(1.65,1.95,121)
# HMoho=np.linspace(20,60,201)
# w1Moho, w2Moho, w3Moho = [0.6, 0.3, 0.1]
# #-------------------------------------------------------------------------
# SynthSeqHKResult=hkSeqSynth(rfSynthSed=rfSynthSed, rfSynthMoho=rfSynthMoho, time=t, rayp=rayp, 
#                              HSed=HSed, KSed=KSed, VpSed=VpSed, w1Sed=w1Sed, w2Sed=w2Sed, w3Sed=w3Sed, 
#                              HMoho=HMoho, KMoho=KMoho, VpMoho=VpMoho, w1Moho=w1Moho, w2Moho=w2Moho, 
#                              w3Moho=w3Moho, g = [75.,10., 15., 2.5], rmneg = None)
# plothkSeqSynth(SynthSeqHKResult, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg', savepath=savepath)

# #%% Test the Synthetic Resonance Filtering
# import numpy as np
# import rffw 
# import os 
# #----------------------------------------------------------
# # Synthetic Parameters
# moddim = 4
# rayp = 0.05826535137278129   #for 60 degree distance
# gaussalp = 1.25
# delay = 10
# n = 2001
# delta = 0.025
# #----------------------------------------------------------
# # Model Parameters
# #----------------------------------------------------------
# depth = np.array([35, 77.5])
# vp = np.array([6.90, 8.045])
# vs = np.array([3.951172240000002, 4.49])
# rho = np.array([2.939564099940003, 3.299])
# rf= rffw.rffw(depth, vp, rho, vs, rayp, gaussalp, delay, n, delta)
# l = len(rf)
# t = np.arange(0, l)
# t = (delta *  t) - delay
# savepath = target_dir1 + '/plots/SynthRF_ResonanceFlt/'
# if not os.path.exists(savepath):  # create data folder if necessary
#     os.mkdir(savepath)
# RemoveResonance= ResonanceFilt(Synthrf=rf, time=t)
# # print(RemoveResonance)
# plotfiltSynthrf(RemoveResonance, savepath, format = 'jpg')
# #%%