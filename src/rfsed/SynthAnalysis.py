"""
.. module:: SynthAnalysis
        :synopsis: Create synthetic reciever function data and analyse the data for hk, sequential hk, 
        and resonance filtering
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""

import numpy as np
from scipy.signal import detrend
from scipy.fft import fft, ifft
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from rfsed.SynRF.FwdRF import SynRF
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from matplotlib.legend_handler import HandlerLine2D as HL
import seaborn
seaborn.set_style("darkgrid", {"axes.facecolor": ".85"})
from scipy import signal
from scipy.signal import argrelextrema, argrelmin, argrelmax
import math
import cmath
#----------------------------------------------------------
def getamp(rfdata, tarray, t):
    """
    Get the amplitude of the receiver function at a specific time

    :param rfdata: receiver function data
    :type rfdata: numpy array
    :param tarray: time array
    :type tarray: numpy array
    :param t: time to get the amplitude
    :type t: float
    :return: amplitude of the receiver function at time t
    """
    amp = rfdata[(np.abs(tarray - t).argmin())]
    return amp
#----------------------------------------------------------
def hkSynth(Synthrf, time, rayp=0.04, H=np.linspace(20,60,201), K=np.linspace(1.65,1.95,121), 
            Vp=6.9, w1=0.6, w2=0.3, w3=0.1, layer = None):
    """
    Calculate the H-K stacking for synthetic receiver function after Zhu and Kanamori (2000)
    Zhu, L., & Kanamori, H. (2000). Moho depth variation in southern california from
    teleseismic receiver functions. Journal of Geophysical Research: Solid Earth, 105 , 
    2969-2980. doi: 10.1029/1999jb900322
    
    :param Synthrf: synthetic receiver function
    :type Synthrf: numpy array
    :param time: time array
    :type time: numpy array
    :param rayp: ray parameter
    :type rayp: float
    :param H: depth array
    :type H: numpy array
    :param K: Vp/Vs array
    :type K: numpy array
    :param Vp: Vp value
    :type Vp: float
    :param w1: weight for Ps
    :param w2: weight for PpPs
    :param w3: weight for PsPs+PpSs
    :type w1, w2, w3: float
    :param layer: layer name 'Moho' or 'Sed'
    :type layer: string
    :return: Dictionary of H-K stacking result for synthetic receiver function
    """
    if Vp is None: Vp = 6.5
    if w1 is None: w1 = 0.6
    if w2 is None: w2 = 0.3
    if w3 is None: w3 = 0.1   
    if layer is None: layer = 'Moho'
    #----------------------------------------------------------
    rp=rayp
    rfdata = Synthrf
    t = time
    stk = np.zeros((len(K)*len(H),3))
    z = 0 
    for i in range(len(K)):
        Ktemp = K[i]
        for j in range(len(H)):
            Htemp = H[j]
            s = 0.0
            trdata = rfdata
            #--------------------------------------
            term1= ((Ktemp/Vp)**2 - rp**2)**0.5
            term2= ((1/Vp)**2 - rp**2)**0.5
            tps = Htemp * (term1 - term2)
            tppps = Htemp * (term1 + term2)
            tpsps = Htemp * 2 * (term1) 
            #--------------------------------------
            ampps = getamp(trdata, t, tps)
            ampppps = getamp(trdata,t, tppps)
            amppsps = getamp(trdata, t, tpsps)
            #--------------------------------------
            stemp = (w1 * ampps) + (w2 * ampppps) - (w3 * amppsps)
            s = stemp + s
            stk[z,:] = [Htemp, Ktemp, s]
            z = z + 1
    bmod = stk[np.argmax(stk[:,2]),:]
    #Show the best model
    plt.tricontourf(stk[:,0], stk[:,1], stk[:,2],60, cmap='jet')
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=plt.plot(bmod[0],bmod[1], 'k+', mew=5, ms=15,\
        label='Best Model %s km %s Vp/Vs'%(bmod[0], bmod[1]))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
        ncol=2, mode="expand", borderaxespad=0, \
        handler_map={p1:HL(numpoints=1)}) 
    plt.ylabel('Vp/Vs')
    plt.xlabel('Depth km')
    plt.xlim(min(H),max(H))
    plt.ylim(min(K),max(K))
    # plt.savefig(savepath , format=format, dpi=250) 
    plt.show()
    plt.close("all")  
    Hbm = bmod[0]
    Kbm = bmod[1]
    Sbm = bmod[2]
    print("Best depth: ", Hbm, "Best Vp/Vs:", Kbm, "Max stack: ", Sbm)
    HKResultSynth = {'SynthRF': Synthrf, 'time':t, 'bestmodel':bmod, 'Vp':Vp, \
                'stackvalues':stk, 'layer':layer, 'rayp':rayp, 'H':H, 'K':K}
    return HKResultSynth
#----------------------------------------------------------
def hkSeqSynth(rfSynthSed, rfSynthMoho, time, rayp=0.04, HSed=np.linspace(0,10,201), KSed=np.linspace(1.65,2.25,201),
               VpSed=3.0, w1Sed=0.6, w2Sed=0.3, w3Sed=0.1, HMoho=np.linspace(20,60,201), KMoho=np.linspace(1.65,1.95,121),
               VpMoho=6.9, w1Moho=0.6, w2Moho=0.3, w3Moho=0.1, g = [75.,10., 15., 2.5], rmneg = None):
    """
    Calculate the Sequential H-K stacking for synthetic receiver function after Yeck et al., 2013
    Yeck, W. L., Sheehan, A. F., & Schulte-Pelkum, V. (2013, 6). Sequential h-k stacking to obtain accurate 
    crustal thicknesses beneath sedimentary basins. Bulletin of the Seismological Society of America, 103 , 2142-2150. doi: 10.1785/0120120290

    :param rfSynthSed: synthetic receiver function for sediment layer (high frequency)
    :param rfSynthMoho: synthetic receiver function for Moho layer (lower frequency)
    :type rfSynthSed, rfSynthMoho: numpy array
    :param time: time array
    :type time: numpy array
    :param rayp: ray parameter
    :type rayp: float
    :param HSed: depth array for sediment layer
    :type HSed: numpy array
    :param KSed: Vp/Vs array for sediment layer
    :type KSed: numpy array
    :param VpSed: Vp value for sediment layer
    :type VpSed: float
    :param w1Sed: weight for Ps
    :param w2Sed: weight for PpPs
    :param w3Sed: weight for PsPs+PpSs
    :type w1Sed, w2Sed, w3Sed: float
    :param HMoho: depth array for Moho layer
    :type HMoho: numpy array
    :param KMoho: Vp/Vs array for Moho layer
    :type KMoho: numpy array
    :param VpMoho: Vp value for Moho layer
    :type VpMoho: float
    :param w1Moho: weight for Ps
    :param w2Moho: weight for PpPs
    :param w3Moho: weight for PsPs+PpSs
    :type w1Moho, w2Moho, w3Moho: float
    :param g: gain for plotting
    :type g: list
    :param rmneg: remove negative values in the stack
    :type rmneg: boolean
    :return: Dictionary of Sequential H-K stacking result
    """
    if g is None: g = [75.,10., 15., 2.5]
    if VpSed is None: VpSed = 3.0
    if w1Sed is None: w1Sed = 0.6
    if w2Sed is None: w2Sed = 0.3
    if w3Sed is None: w3Sed = 0.1  
    if rmneg is None: rmneg = False
    rfdata = rfSynthMoho
    rp=rayp
    t = time
    #----------------------------------------------------------
    # Stack for Sediment Layer
    #----------------------------------------------------------
    stkSed = np.zeros((len(KSed)*len(HSed),3))
    z = 0 
    for i in range(len(KSed)):
        Ktemp = KSed[i]
        for j in range(len(HSed)):
            Htemp = HSed[j]
            s = 0.0
            trdata = rfSynthSed
            rp = rp
            #--------------------------------------
            term1= ((Ktemp/VpSed)**2 - rp**2)**0.5
            term2= ((1/VpSed)**2 - rp**2)**0.5
            tpsSed = Htemp * (term1 - term2)
            tpppsSed = Htemp * (term1 + term2)
            tpspsSed = Htemp * 2 * (term1) 
            #--------------------------------------
            ampps = getamp(trdata, t, tpsSed)
            ampppps = getamp(trdata,t, tpppsSed)
            amppsps = getamp(trdata, t, tpspsSed)
            #--------------------------------------
            stemp = (w1Sed * ampps) + (w2Sed * ampppps) - (w3Sed * amppsps)
            sSed = stemp + s
            stkSed[z,:] = [Htemp, Ktemp, sSed]
            z = z + 1
    bmodSed = stkSed[np.argmax(stkSed[:,2]),:]
    #----------------------------------------------------------
    # Plot Sediment Stacking results
    plt.tricontourf(stkSed[:,0], stkSed[:,1], stkSed[:,2],60, cmap='jet')
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=plt.plot(bmodSed[0],bmodSed[1], 'k+', mew=5, ms=15,\
        label='Best Model %s km %s Vp/Vs'%(bmodSed[0], bmodSed[1]))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
        ncol=2, mode="expand", borderaxespad=0, \
        handler_map={p1:HL(numpoints=1)}) 
    plt.ylabel('Vp/Vs')
    plt.xlabel('Depth km')
    plt.xlim(min(HSed),max(HSed))
    plt.ylim(min(KSed),max(KSed))
    # plt.savefig(savepath , format=format, dpi=250) 
    plt.show()
    plt.close("all")    
    #----------------------------------------------------------
    # Add Ps, PpPs, PsPs + PbSs timing to the data header
    #----------------------------------------------------------    
    H1 = bmodSed[0]
    K1 = bmodSed[1]
    SSed = bmodSed[2]
    print("Best Sediment depth: ", H1, "Best Sediment Vp/Vs:", K1, "Max stack: ", SSed)
    #----------------------------------------------------------
    # Stack for Moho Layer
    #----------------------------------------------------------
    stkMoho = np.zeros((len(KMoho)*len(HMoho),3))
    z = 0 
    for i in range(len(KMoho)):
        Ktemp = KMoho[i]
        for j in range(len(HMoho)):
            Htemp = HMoho[j]
            s = 0.0
            trdata = rfSynthMoho
            rp = rp
            #--------------------------------------
            term1 = ((K1/VpSed)**2-(rp)**2)**0.5
            term2 = ((1/VpSed)**2-(rp)**2)**0.5
            term3 = ((Ktemp/VpMoho)**2-(rp)**2)**0.5
            term4 = ((1/VpMoho)**2-(rp)**2)**0.5
            #--------------------------------------
            tpsMoho = (H1 * (term1 - term2)) + (Htemp * (term3 - term4))
            tpppsMoho = (H1 * (term1 + term2)) + (Htemp * (term3 + term4))
            tpspsMoho = (2 * H1 * term1) + (2 * Htemp * term3)
            #--------------------------------------
            ampps = getamp(trdata, t, tpsMoho) 
            ampppps = getamp(trdata, t, tpppsMoho)  
            amppsps = getamp(trdata, t, tpspsMoho)
            stemp = (w1Moho * ampps) + (w2Moho * ampppps) - (w3Moho * amppsps)
            sMoho = stemp + s
                #--------------------------------------
            stkMoho[z,:] = [Htemp, Ktemp, sMoho]
            z = z + 1
    bmodMoho = stkMoho[np.argmax(stkMoho[:,2]),:]
    #----------------------------------------------------------
    # plot Moho Stacking results
    for i in range(len(stkMoho)):
        if stkMoho[i,2] <= 0:
            stkMoho[i,2] = 0
    plt.tricontourf(stkMoho[:,0], stkMoho[:,1], stkMoho[:,2],60, cmap='jet')
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=plt.plot(bmodMoho[0],bmodMoho[1], 'k+', mew=5, ms=15,\
        label='Best Model %s km %s Vp/Vs'%(bmodMoho[0], bmodMoho[1]))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
        ncol=2, mode="expand", borderaxespad=0, \
        handler_map={p1:HL(numpoints=1)}) 
    plt.ylabel('Vp/Vs')
    plt.xlabel('Depth km')
    plt.xlim(min(HMoho),max(HMoho))
    plt.ylim(min(KMoho),max(KMoho))
    # plt.savefig(savepath , format=format, dpi=250) 
    plt.show()
    plt.close("all")    
    H2 = bmodMoho[0]
    K2 = bmodMoho[1] 
    print("Best Moho depth: ", H2, "Best Moho Vp/Vs:", K2, "Max stack: ", SSed)   
    SynthSeqHKResult = {'SynthRFSed': rfSynthSed, 'SynthRFMoho': rfSynthMoho, 'time':t, \
                          'bestmodelSed':bmodSed, 'stackvalueSed':stkSed, 'bestmodelMoho':bmodMoho, 'stackvalueMoho':stkMoho,
                            'RayP':rayp, 'VpSed':VpSed, 'VpMoho':VpMoho, 'HSed':HSed, 'KSed':KSed, 'HMoho':HMoho, 'KMoho':KMoho}
    return SynthSeqHKResult
#----------------------------------------------------------
def ResonanceFilt(Synthrf, time):
    """
    Resonance Filter for synthetic receiver function data (Removes sediment reverberation effect)

    :param Synthrf: synthetic receiver function data
    :type Synthrf: numpy array
    :param time: time array
    :type time: numpy array

    :return: A dictionary of input synthetic receiver function, Filtered receiver function, 
            autocorrelation of the data, resonance filter, time lag (2 way travel time of the sediment reverbration), 
            and the strength of the sediment reverbration (r)
    """
    t = time
    delta=t[1]-t[0]
    dt=delta
    #---------------------------------------------
    # Build the frequency vector
    #---------------------------------------------
    n = len(Synthrf)
    fmax = 1 / (2.0 * dt)
    df = fmax / (n / 2)
    f = np.hstack((df * np.arange(0, n//2+1), df * np.arange(-n//2 + 1, 0)))
    nf = n // 2 + 1
    dw = 2.0 * np.pi * df
    w = dw * np.hstack((np.arange(0, n//2+1), np.arange(-n//2 + 1, 0)))
    filtered_rf = np.zeros_like(Synthrf)
    D = Synthrf
    D = D - np.mean(D)
    D = detrend(D)
    #-----------------------------------------------
    # Calculate the autocorrelation
    ac = np.correlate(D, D, mode='full')
    ac = ac / np.max(ac)
    ac = ac[n-1:2*n-1]

    rac = -ac
    locs, _ = find_peaks(rac)
    if len(locs) == 0:
        r0 = 0
        tlag = 0
        print('No reverberation detected')
    else:
        r0 = np.abs(rac[locs[0]])
        tlag = locs[0] * dt

    resonanceflt = (1 + r0 * np.exp(-1j * w * tlag))
    filtered_rf = np.real(ifft(fft(D) * resonanceflt))
    # filtered_rf = filtered_rf / np.max(filtered_rf)
    #------------------------------------------------------------------------------
    FltResults={'rf':Synthrf, 'filteredrf': filtered_rf, 'resonancefilter':resonanceflt, 'time':time, 'autoc':ac, 'r':r0, 'tlag':tlag, 'delta':delta}
    return FltResults

def plothkSynth(HKResultSynth, savepath, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg'): 
    """
    Plot the H-K stacking result for synthetic receiver function

    :param HKResultSynth: Dictionary containing H-K stacking result from function hkSynth
    :type HKResultSynth: dictionary
    :param g: gain for plotting
    :type g: list
    :param rmneg: remove negative values in the stack
    :type rmneg: boolean
    :param format: format for saving the figure
    :type format: string
    :param savepath: path to save the figure
    :type savepath: string
    """
    if rmneg is None: rmneg = False
    if format is None: format = 'pdf'
    if g is None: g = [75.,10., 15., 2.5]
    staname = 'SynthRF'
    bmod = HKResultSynth['bestmodel']
    Hbm = bmod[0]
    Kbm = bmod[1]
    Sbm = bmod[2]
    layer = HKResultSynth['layer']
    rfstream=HKResultSynth['SynthRF']
    rp = HKResultSynth['rayp']
    t=HKResultSynth['time']
    Vp=HKResultSynth['Vp']
    term1= ((Kbm/Vp)**2 - rp**2)**0.5
    term2= ((1/Vp)**2 - rp**2)**0.5
    tps = Hbm * (term1 - term2)
    tppps = Hbm * (term1 + term2)
    tpsps = Hbm * 2 * (term1)
    #----------------------------------------------------------
    # plot H-K results
    #----------------------------------------------------------
    if HKResultSynth['layer'] == 'Moho':
        mintime=-5
        maxtime=30
        span=0.5
        timestep=5
    elif HKResultSynth['layer'] == 'Sed':
        mintime=-1
        maxtime=10
        span=0.1
        timestep=1
    stk = HKResultSynth['stackvalues']
    H = HKResultSynth['H']
    K = HKResultSynth['K']
    plt.figure(figsize=(18, 6)) 
    gs = gridspec.GridSpec(1, 3, width_ratios=[6,4,4])
    #-----Remove the stack values that is lower than zero------
    if rmneg == True:
        for i in range(len(stk)):
            if stk[i,2] <= 0:
                stk[i,2] = 0
    #----------------------------------------------------------
    ax1 = plt.subplot(gs[0, 0])
    plt.tricontourf(stk[:,0], stk[:,1], stk[:,2],50, cmap ="jet")
    #----------------------------------------------------------
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=ax1.plot(bmod[0],bmod[1], 'k+', mew=2, ms=10, \
    label='%s Best Model %.2f km %.2f Vp/Vs'%(staname,bmod[0],bmod[1]+0.001))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0, handler_map={p1:HL(numpoints=1)})
    ax1.set_ylabel('Vp/Vs')
    ax1.set_xlabel('Depth km')
    ax1.set_xlim(min(H), max(H))
    ax1.set_ylim(min(K), max(K))
    #----------------------------------------------------------
    # Plot the receiver function by Back Azimuth with the estimated times
    #----------------------------------------------------------
    ax2 = plt.subplot(gs[0, 1]) 
    gcarcarr = np.zeros(len(rfstream))
    baz = 180
    rfdata = ((rfstream) * g[0]) + baz
    t = t
    major_ticks_x = np.arange(-10, 41, timestep)                                              
    minor_ticks_x = np.arange(-10, 41, 1) 
    major_ticks_y = np.arange(0, 361, 45)                                              
    minor_ticks_y = np.arange(-30, 390, 5)
    ax2.set_xticks(major_ticks_x)                                                       
    ax2.set_xticks(minor_ticks_x, minor=True)                                           
    ax2.set_yticks(major_ticks_y)                                                       
    ax2.set_yticks(minor_ticks_y, minor=True)                                          
    ax2.set_xlim(mintime, maxtime)
    ax2.set_ylim(-29, 390)
    ax2.plot(t, rfdata, "k-", lw=0.5)
    ax2.fill_between(t, baz, rfdata, where=rfdata > baz, facecolor='red', alpha = 0.25)
    ax2.fill_between(t, baz, rfdata, where=rfdata <= baz, facecolor='blue', alpha = 0.25)
    mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    #----------------------------------------------------------
    # Plot H-K sequential formula times
    #----------------------------------------------------------
    ampgain = g[1]
    ax2.plot([tps, tps], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    ax2.plot([tppps, tppps], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    ax2.plot([tpsps, tpsps], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    #---------------------------------------------------------
    # Highlight specifiec time range on RF plot
    #----------------------------------------------------------
    ax2.axvspan(tps-span, tps+span, alpha=0.5, color='lightgreen')
    ax2.axvspan(tppps-span, tppps+span, alpha=0.5, color='lightgreen')
    ax2.axvspan(tpsps-span, tpsps+span, alpha=0.5, color='lightgreen')
    ax2.set_xlabel("Time (sec)")
    ax2.set_ylabel("Back-Azimuth (deg)")    
    #----------------------------------------------------------
    filename = "%s/%s_%s.%s"%(savepath, staname, str(layer), format)
    plt.savefig(filename , format=format, transparent=False,\
        dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close('all')
#----------------------------------------------------------
def plothkSeqSynth(SynthSeqHKResult, savepath, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg'): 
    """
    Plot the Sequential H-K stacking result for synthetic receiver function

    :param SynthSeqHKResult: Dictionary containing Sequential H-K stacking result from function hkSeqSynth
    :type SynthSeqHKResult: dictionary
    :param g: gain for plotting
    :type g: list
    :param rmneg: remove negative values in the stack
    :type rmneg: boolean
    :param format: format for saving the figure
    :type format: string
    :param savepath: path to save the figure
    :type savepath: string
    """
    if rmneg is None: rmneg = False
    if format is None: format = 'pdf'
    if g is None: g = [75.,10., 15., 2.5]
    staname='SynthRF_SeqHK'
    VpSed=SynthSeqHKResult['VpSed']
    HSed=SynthSeqHKResult['HSed']
    KSed=SynthSeqHKResult['KSed']
    bmodSed = SynthSeqHKResult['bestmodelSed']
    H1 = bmodSed[0]
    K1 = bmodSed[1]
    stkSed=SynthSeqHKResult['stackvalueSed']
    
    rfstream=SynthSeqHKResult['SynthRFSed']
    rp = SynthSeqHKResult['RayP']
    t=SynthSeqHKResult['time']
    term1= ((K1/VpSed)**2 - rp**2)**0.5
    term2= ((1/VpSed)**2 - rp**2)**0.5
    tpsSed = H1 * (term1 - term2)
    tpppsSed = H1 * (term1 + term2)
    tpspsSed = H1 * 2 * (term1)
    #----------------------------------------------------------
    # plot H-K results
    #----------------------------------------------------------
    plt.figure(figsize=(18, 6)) 
    gs = gridspec.GridSpec(1, 3, width_ratios=[6,4,4])
    #----------------------------------------------------------
    # plot misfit
    #----------------------------------------------------------
    #-----Remove the stack values that is lower than zero------
    if rmneg == True:
        for i in range(len(stkSed)):
            if stkSed[i,2] <= 0:
                stkSed[i,2] = 0
    #----------------------------------------------------------
    ax1 = plt.subplot(gs[0, 0])
    plt.tricontourf(stkSed[:,0], stkSed[:,1], stkSed[:,2],50, cmap ="jet")
    #----------------------------------------------------------
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=ax1.plot(bmodSed[0],bmodSed[1], 'k+', mew=2, ms=10, \
    label='%s Best Model %.2f km %.2f Vp/Vs'%(staname,bmodSed[0],bmodSed[1]+0.001))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0, handler_map={p1:HL(numpoints=1)})
    ax1.set_ylabel('Vp/Vs')
    ax1.set_xlabel('Depth km')
    ax1.set_xlim(min(HSed), max(HSed))
    ax1.set_ylim(min(KSed), max(KSed))
    #----------------------------------------------------------
    # Plot the receiver function witht the estimated times
    #----------------------------------------------------------
    ax2 = plt.subplot(gs[0, 1]) 
    baz=180
    rfdata = ((SynthSeqHKResult['SynthRFSed']) * (g[0]/2)) + baz
    major_ticks_x = np.arange(-10, 41, 1)                                              
    minor_ticks_x = np.arange(-10, 41, 1) 
    major_ticks_y = np.arange(0, 361, 45)                                              
    minor_ticks_y = np.arange(-30, 390, 5)
    ax2.set_xticks(major_ticks_x)                                                       
    ax2.set_xticks(minor_ticks_x, minor=True)                                           
    ax2.set_yticks(major_ticks_y)                                                       
    ax2.set_yticks(minor_ticks_y, minor=True)                                          
    ax2.set_xlim(-1, 10)
    ax2.set_ylim(-29, 390)
    ax2.set_xlabel("Time (sec)")
    ax2.set_ylabel("Back-Azimuth (deg)")  
    ax2.plot(t, rfdata, "k-", lw=0.5)
    ax2.fill_between(t, baz, rfdata, where=rfdata > baz, facecolor='red', alpha = 0.25)
    ax2.fill_between(t, baz, rfdata, where=rfdata <= baz, facecolor='blue', alpha = 0.25)
    mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    ampgain = g[1]
    ax2.plot([tpsSed, tpsSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5, label="Moho")
    ax2.plot([tpppsSed, tpppsSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    ax2.plot([tpspsSed, tpspsSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    filename = "%s/%s_%s.%s"%(savepath, staname, "Sediment", format)
    plt.savefig(filename , format=format, transparent=False,\
        dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close('all')
    #----------------------------------------------------------
    # Plot Moho Stack
    #----------------------------------------------------------
    VpMoho=SynthSeqHKResult['VpMoho']
    HMoho=SynthSeqHKResult['HMoho']
    KMoho=SynthSeqHKResult['KMoho']
    bmodMoho = SynthSeqHKResult['bestmodelMoho']  
    H2 = bmodMoho[0]
    K2 = bmodMoho[1]
    stkMoho=SynthSeqHKResult['stackvalueMoho']
    rfstream=SynthSeqHKResult['SynthRFMoho']
    rp = rp
    t=SynthSeqHKResult['time']
    term1 = ((K1/VpSed)**2-(rp)**2)**0.5    
    term2 = ((1/VpSed)**2-(rp)**2)**0.5
    term3 = ((K2/VpMoho)**2-(rp)**2)**0.5    
    term4 = ((1/VpMoho)**2-(rp)**2)**0.5
    #--------------------------------------
    tpsMoho = (H1 * (term1 - term2)) + (H2 * (term3 - term4))
    tpppsMoho = (H1 * (term1 + term2)) + (H2 * (term3 + term4))
    tpspsMoho = (2 * H1 * term1) + (2 * H2 * term3)
    #----------------------------------------------------------
    plt.figure(figsize=(18, 6)) 
    gs = gridspec.GridSpec(1, 3, width_ratios=[6,4,4])
    #-----Remove the stack values that is lower than zero------
    if rmneg == True:
        for i in range(len(stkMoho)):
            if stkMoho[i,2] <= 0:
                stkMoho[i,2] = 0
    #----------------------------------------------------------
    ax1 = plt.subplot(gs[0, 0])
    plt.tricontourf(stkMoho[:,0], stkMoho[:,1], stkMoho[:,2],50, cmap ="jet")
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    cb.ax.set_xlabel('S', rotation=0)
    p1,=ax1.plot(bmodMoho[0],bmodMoho[1], 'k+', mew=2, ms=10, \
    label='%s Best Model %.2f km %.2f Vp/Vs'%(staname,bmodMoho[0],bmodMoho[1]+0.001))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
        ncol=2, mode="expand", borderaxespad=0, \
        handler_map={p1:HL(numpoints=1)})
    ax1.set_ylabel('Vp/Vs')
    ax1.set_xlabel('Depth km')
    ax1.set_xlim(min(HMoho), max(HMoho))
    ax1.set_ylim(min(KMoho), max(KMoho))
    #----------------------------------------------------------
    # Plot the receiver function witht the estimated times
    #----------------------------------------------------------
    ax2 = plt.subplot(gs[0, 1]) 
    baz=180
    rfdata = ((SynthSeqHKResult['SynthRFMoho']) * g[0]) + baz
    major_ticks_x = np.arange(-10, 41, 5)                                           
    minor_ticks_x = np.arange(-10, 41, 1) 
    major_ticks_y = np.arange(0, 361, 45)                                              
    minor_ticks_y = np.arange(-30, 390, 5)
    ax2.set_xticks(major_ticks_x)                                                       
    ax2.set_xticks(minor_ticks_x, minor=True)                                           
    ax2.set_yticks(major_ticks_y)                                                       
    ax2.set_yticks(minor_ticks_y, minor=True)                                  
    ax2.set_xlim(-5, 30)
    ax2.set_ylim(-29, 390)
    ax2.plot(t, rfdata, "k-", lw=0.5)
    ax2.fill_between(t, baz, rfdata, where=rfdata > baz, 
                        facecolor='red', alpha = 0.25)
    ax2.fill_between(t, baz, rfdata, where=rfdata <= baz, 
                        facecolor='blue', alpha = 0.25)
    mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    #----------------------------------------------------------
    # Plot H-K sequential formula times
    #----------------------------------------------------------
    ampgain = g[1]
    ax2.plot([tpsMoho, tpsMoho], [baz+ampgain, baz-ampgain], 'g-', lw=1.5, label="Moho")
    ax2.plot([tpppsMoho, tpppsMoho], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    ax2.plot([tpspsMoho, tpspsMoho], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    filename = "%s/%s_%s.%s"%(savepath, staname, 'Moho', format)
    plt.savefig(filename , format=format, transparent=False,\
    dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close('all')
#----------------------------------------------------------
def plotfiltSynthrf(FilteredRF, savepath, format = 'jpg'):
    """
    Plot the filtered receiver function, autocorrelation, and the resonance filter

    :param FilteredRF: Dictionary of results from the resonance filtering method (function: Resonance_Filt)
    :type FilteredRF: dictionary
    :param savepath: path to save the plots
    :type savepath: str
    :param format: format of the plot (default: jpg)
    :type format: str
    """
    rf=FilteredRF['rf']
    filtered_rf=FilteredRF['filteredrf']
    resonanceflt=FilteredRF['resonancefilter']
    time=FilteredRF['time']
    autoc=FilteredRF['autoc']
    suff='SyntheticRF'
    # Plot Autocorrelation
    plt.plot(time,autoc)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Autocorrelation of %s'%suff)
    filename = "%s/%s.%s"%(savepath, 'Autocorrelation', format)
    plt.savefig(filename , format=format, transparent=False,\
            dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
    # Plot Resonance Filter
    plt.plot(time, resonanceflt)
    plt.xlabel('Time (s)')
    # plt.ylabel('Amplitude')
    plt.title('Resonance Filter of %s'%suff)
    filename = "%s/%s.%s"%(savepath, "Resonance_Filter", format)
    plt.savefig(filename , format=format, transparent=False,\
            dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
    # Plot RF and Filtered RF
    plt.plot(time, rf, label = 'RF', color = 'black')
    plt.plot(time, filtered_rf, label = 'Filtered RF', color = 'red')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('RF and Filtered RF of %s'%suff)
    plt.legend()
    filename = "%s/%s.%s"%(savepath, "RF_FilteredRF", format)
    plt.savefig(filename , format=format, transparent=False,\
            dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
#----------------------------------------------------------