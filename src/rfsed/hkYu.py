"""
.. module:: hkYu
        :synopsis: Modified H-K stacking method of Yu et al. (2015) for receiver functions filtered with the resonance filter
        Yu, Y., Song, J., Liu, K. H., & Gao, S. S. (2015). Determining crustal structure
        beneath seismic stations overlying a low-velocity sedimentary layer using receiver
        functions. Journal of Geophysical Research: Solid Earth, 120 , 3208-3218. doi:
        10.1002/2014JB011610
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""

import numpy as np
from scipy import signal
from scipy.signal import argrelextrema, argrelmin, argrelmax
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from scipy.fftpack import fft, ifft
from matplotlib.legend_handler import HandlerLine2D as HL
import seaborn
seaborn.set_style("darkgrid", {"axes.facecolor": ".85"})

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
def hkYu(FltResult, rayp=0.04, HSubSed=np.linspace(20,60,201), KSubSed=np.linspace(1.65,1.95,121), 
               HSed=np.linspace(0,10,201), KSed=np.linspace(1.65,2.25,201), VpMoho=6.9, VpSed= 2.5,  
               w1SubSed=0.6, w2SubSed=0.3, w3SubSed=0.1, w1Sed=0.6, w2Sed=0.3, w3Sed=0.1):
    """
    Modified H-K stacking method of Yu et al. (2015)

    :param FltResult: Dictionary of results from the resonance filtering method (function: Resonance_Filt)
    :type FltResult: dict
    :param rayp: Ray parameter
    :type rayp: float
    :param HSubSed: Subsediment layer thickness array
    :type HSubSed: numpy array
    :param KSubSed: Subsediment layer Vp/Vs array
    :type KSubSed: numpy array
    :param HSed: Sediment layer thickness array
    :type HSed: numpy array
    :param KSed: Sediment layer Vp/Vs array
    :type KSed: numpy array
    :param VpMoho: Moho Vp
    :type VpMoho: float
    :param VpSed: Sediment Vp
    :type VpSed: float
    :param w1SubSed: Weight for the Subsediment Ps arrival at adjusted arrival time
    :type w1SubSed: float
    :param w2SubSed: Weight for the Subsediment Ppps arrival at adjusted arrival time
    :type w2SubSed: float
    :param w3SubSed: Weight for the Subsediment Psps+PpSs arrival at adjusted arrival time
    :type w3SubSed: float
    :param w1Sed: Weight for the Ps of the sediment layer
    :type w1Sed: float
    :param w2Sed: Weight for the Ppps of the Subsediment layer
    :type w2Sed: float
    :param w3Sed: Weight for the Psps+PpSs of the Subsediment layer
    :type w3Sed: float
    """
    if VpMoho is None: VpMoho = 6.5
    if VpSed is None: VpSed = 2.5
    if w1SubSed is None: w1 = 0.6
    if w2SubSed is None: w2 = 0.3
    if w3SubSed is None: w3 = 0.1   
    #----------------------------------------------------------
    rp=rayp
    Fltrf=FltResult['filteredrf']
    rfdata = Fltrf
    t = FltResult['time']
    tlag=FltResult['tlag']
    VsSed=0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
    SedH=tlag*VsSed/2
    Pdelay=SedH/VsSed - SedH/VpSed
    #----------------------------------------------------------
    # Stack for SubSediment Layer
    stkSubSed = np.zeros((len(KSubSed)*len(HSubSed),3))
    z = 0 
    for i in range(len(KSubSed)):
        Ktemp = KSubSed[i]
        for j in range(len(HSubSed)):
            Htemp = HSubSed[j]
            s = 0.0
            trdata = rfdata
            #--------------------------------------
            term1= ((Ktemp/VpMoho)**2 - rp**2)**0.5
            term2= ((1/VpMoho)**2 - rp**2)**0.5
            tps = Htemp * (term1 - term2)
            tppps = Htemp * (term1 + term2)
            tpsps = Htemp * 2 * (term1)
            #--------------------------------------
            # Adjusted Arrival Times for subsediment layer stacking
            tpsSubSed=tps+Pdelay
            tpppsSubSed=tppps+tlag - Pdelay
            tpspsSubSed=tpsps+tlag
            #--------------------------------------
            ampps = getamp(trdata, t, tpsSubSed)
            ampppps = getamp(trdata,t, tpppsSubSed)
            amppsps = getamp(trdata, t, tpspsSubSed)
            #--------------------------------------
            stemp = (w1SubSed * ampps) + (w2SubSed * ampppps) - (w3SubSed * amppsps)
            s = stemp + s
            stkSubSed[z,:] = [Htemp, Ktemp, s]
            z = z + 1
    bmodSubSed = stkSubSed[np.argmax(stkSubSed[:,2]),:]
    #Show the best model
    plt.tricontourf(stkSubSed[:,0], stkSubSed[:,1], stkSubSed[:,2],60, cmap='jet')
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=plt.plot(bmodSubSed[0],bmodSubSed[1], 'k+', mew=5, ms=15,\
        label='Best Model %s km %s Vp/Vs'%(bmodSubSed[0], bmodSubSed[1]))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
        ncol=2, mode="expand", borderaxespad=0, \
        handler_map={p1:HL(numpoints=1)}) 
    plt.ylabel('Vp/Vs')
    plt.xlabel('Depth km')
    plt.xlim(min(HSubSed),max(HSubSed))
    plt.ylim(min(KSubSed),max(KSubSed))
    # plt.savefig(savepath , format=format, dpi=250) 
    plt.show()
    plt.close("all")  
    HSubSedbm = bmodSubSed[0]
    KSubSedbm = bmodSubSed[1]
    SSubSedbm = bmodSubSed[2]
    print("Best Subsediment thickness: ", HSubSedbm, "Best Subsediment Vp/Vs:", KSubSedbm, "Max stack: ", SSubSedbm)
    #----------------------------------------------------------
    # Stack for Sediment Layer
    stkSed = np.zeros((len(KSed)*len(HSed),3))
    z = 0 
    for i in range(len(KSed)):
        Ktemp = KSed[i]
        for j in range(len(HSed)):
            Htemp = HSed[j]
            s = 0.0
            trdata = rfdata
            #--------------------------------------
            # Calculate tps for Sediment phase
            term3= ((Ktemp/VpSed)**2 - rp**2)**0.5
            term4= ((1/VpSed)**2 - rp**2)**0.5
            tpsSed = Htemp * (term3 - term4)
            #--------------------------------------
            # tpps and tpsps are fixed based on the best subsediment model
            term5= ((KSubSedbm /VpMoho)**2 - rp**2)**0.5
            term6= ((1/VpMoho)**2 - rp**2)**0.5
            tpppsSubSed = HSubSedbm * (term5 + term6)
            tpspsSubSed = HSubSedbm * 2 * (term5)  
            #--------------------------------------
            ampps = getamp(trdata, t, tpsSed)
            ampppps = getamp(trdata,t, tpppsSubSed)
            amppsps = getamp(trdata, t, tpspsSubSed)
            #--------------------------------------
            stemp = (w1SubSed * ampps) + (w2SubSed * ampppps) - (w3SubSed * amppsps)
            s = stemp + s
            stkSed[z,:] = [Htemp, Ktemp, s]
            z = z + 1
    bmodSed = stkSed[np.argmax(stkSed[:,2]),:]
    #Show the best model
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
    HSedbm = bmodSed[0]
    KSedbm = bmodSed[1]
    SSedbm = bmodSed[2]
    print("Best Sediment thickness: ", HSedbm, " km", "Best Sediment Vp/Vs:", KSedbm, "Max stack: ", SSedbm)
    print("Best Subsediment thickness: ", HSubSedbm, " km", "Best Subsediment Vp/Vs:", KSubSedbm, "Max stack: ", SSubSedbm)
    HMohobm=HSubSedbm + HSedbm
    print("Best Moho thickness: ", HMohobm, " km")
    HKYuResult = {'FltRF': Fltrf, 'time':t, 'tlag':tlag, 'Pdelay':Pdelay, 'rayp':rayp, 'BestMohoDepth': HMohobm, 
                  'bestmodelSubSed':bmodSubSed, 'VpMoho':VpMoho, 'VpSed':VpSed, 'stackvaluesSubSed':stkSubSed, 
                  'bestmodelSed':bmodSed, 'stackvaluesSed':stkSed,'HMoho':HSubSed, 'KMoho':KSubSed,'HSed':HSed, 
                  'KSed':KSed,'SubSedThick':HSubSedbm, 'SubSedVpVs':KSubSedbm, 'SedThick':HSedbm, 'SedVpVs':KSedbm}
    return HKYuResult
#----------------------------------------------------------
def plothkYu(hkYuResult, savepath, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg'): 
    """
    Plot the results of the Modified H-K stacking method of Yu et al. (2015)

    :param hkYuResult: Dictionary of results from the Modified H-K stacking method (function: hkYu)
    :type hkYuResult: dict
    :param savepath: path to save the plots
    :type savepath: str
    :param g: gain values for the plot (default: [75.,10., 15., 2.5])
    :type g: list
    :param rmneg: remove the stack values that is lower than zero (default: False)
    :type rmneg: bool
    :param format: format of the plot (default: jpg)
    :type format: str   
    """
    if rmneg is None: rmneg = False
    if format is None: format = 'pdf'
    if g is None: g = [75.,10., 15., 2.5]
    staname='FilteredRF_HKYu'
    Pdelay=hkYuResult['Pdelay']
    tlag=hkYuResult['tlag']
    VpMoho=hkYuResult['VpMoho']
    VpSed=hkYuResult['VpSed']
    HSed=hkYuResult['HSed']
    KSed=hkYuResult['KSed']
    HMoho=hkYuResult['HMoho']
    KMoho=hkYuResult['KMoho']
    bmodSed = hkYuResult['bestmodelSed']
    H1 = bmodSed[0]
    K1 = bmodSed[1]
    stkSed=hkYuResult['stackvaluesSed']
    bmodSubSed = hkYuResult['bestmodelSubSed']
    H2 = bmodSubSed[0]
    K2 = bmodSubSed[1]
    stkSubSed=hkYuResult['stackvaluesSubSed']
    #----------------------------------------------------------
    # Plot for Sediment Layer
    rfdata=hkYuResult['FltRF']
    rp = hkYuResult['rayp']
    t=hkYuResult['time']
    #--------------------------------------
    # Calculate tps for Sediment phase
    term3= ((K1/VpSed)**2 - rp**2)**0.5
    term4= ((1/VpSed)**2 - rp**2)**0.5
    tpsSed = H1 * (term3 - term4)
    #--------------------------------------
    # tpps and tpsps are fixed based on the best subsediment model
    term5= ((K2 /VpMoho)**2 - rp**2)**0.5
    term6= ((1/VpMoho)**2 - rp**2)**0.5
    tpppsSubSed = H2 * (term5 + term6)
    tpspsSubSed = H2 * 2 * (term5)  
    #----------------------------------------------------------
    # plot H-K results
    #----------------------------------------------------------
    plt.figure(figsize=(18, 6)) 
    gs = gridspec.GridSpec(1, 3, width_ratios=[6,4,4])
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
    label='Best Model Sediment %.2f km %.2f Vp/Vs'%(bmodSed[0],bmodSed[1]+0.001))
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
    rfdata = ((rfdata) * (g[0]/2)) + baz
    major_ticks_x = np.arange(-10, 41, 5)                                              
    minor_ticks_x = np.arange(-10, 41, 5) 
    major_ticks_y = np.arange(0, 361, 45)                                              
    minor_ticks_y = np.arange(-30, 390, 5)
    ax2.set_xticks(major_ticks_x)                                                       
    ax2.set_xticks(minor_ticks_x, minor=True)                                           
    ax2.set_yticks(major_ticks_y)                                                       
    ax2.set_yticks(minor_ticks_y, minor=True)                                          
    ax2.set_xlim(-5, 30)
    ax2.set_ylim(-29, 390)
    ax2.set_xlabel("Time (sec)")
    ax2.set_ylabel("Back-Azimuth (deg)")  
    ax2.plot(t, rfdata, "k-", lw=0.5)
    ax2.fill_between(t, baz, rfdata, where=rfdata > baz, facecolor='red', alpha = 0.25)
    ax2.fill_between(t, baz, rfdata, where=rfdata <= baz, facecolor='blue', alpha = 0.25)
    mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    ampgain = g[1]
    ax2.plot([tpsSed, tpsSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5, label="Moho")
    ax2.plot([tpppsSubSed, tpppsSubSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    ax2.plot([tpspsSubSed, tpspsSubSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    filename = "%s/%s_%s.%s"%(savepath, staname, "Sediment", format)
    plt.savefig(filename , format=format, transparent=False,\
        dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close('all')
    #----------------------------------------------------------
    # Plot SubSediment Stack
    #----------------------------------------------------------
    rfdata=hkYuResult['FltRF']
    bmodSubSed = hkYuResult['bestmodelSubSed']
    H2 = bmodSubSed[0]
    K2 = bmodSubSed[1]
    stkSubSed=hkYuResult['stackvaluesSubSed']
    term1= ((K2/VpMoho)**2 - rp**2)**0.5
    term2= ((1/VpMoho)**2 - rp**2)**0.5
    tps = H2 * (term1 - term2)
    tppps = H2 * (term1 + term2)
    tpsps = H2 * 2 * (term1)
    #--------------------------------------
    # Adjusted Arrival Times for subsediment layer stacking
    tpsSubSed=tps+Pdelay
    tpppsSubSed=tppps+tlag - Pdelay
    tpspsSubSed=tpsps+tlag
    #----------------------------------------------------------
    plt.figure(figsize=(18, 6)) 
    gs = gridspec.GridSpec(1, 3, width_ratios=[6,4,4])
    #-----Remove the stack values that is lower than zero------
    if rmneg == True:
        for i in range(len(stkSubSed)):
            if stkSubSed[i,2] <= 0:
                stkSubSed[i,2] = 0
    #----------------------------------------------------------
    ax1 = plt.subplot(gs[0, 0])
    plt.tricontourf(stkSubSed[:,0], stkSubSed[:,1], stkSubSed[:,2],50, cmap ="jet")
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    cb.ax.set_xlabel('S', rotation=0)
    p1,=ax1.plot(bmodSubSed[0],bmodSubSed[1], 'k+', mew=2, ms=10, \
    label='Best SubSediment Model %.2f km %.2f Vp/Vs'%(bmodSubSed[0],bmodSubSed[1]+0.001))
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
    rfdata = ((rfdata) * g[0]) + baz
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
    ax2.plot([tpsSubSed, tpsSubSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5, label="Moho")
    ax2.plot([tpppsSubSed, tpppsSubSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    ax2.plot([tpspsSubSed, tpspsSubSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    filename = "%s/%s_%s.%s"%(savepath, staname, 'SubSediment', format)
    plt.savefig(filename , format=format, transparent=False,\
    dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close('all')
#----------------------------------------------------------
