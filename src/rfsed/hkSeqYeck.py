# Copyright (c) 2023, Stephen Akinremi

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
.. module:: hkSeqYeck
        :synopsis: Sequential H-K stacking method for receiver function 
        analysis after Yeck et al., 2013: 
        Yeck, W. L., Sheehan, A. F., & Schulte-Pelkum, V. (2013, 6).
        Sequential h-k stacking to obtain accurate crustal thicknesses beneath
        sedimentary basins. Bulletin of the Seismological Society of America,
        103, 2142-2150. doi: 10.1785/0120120290
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> 
                  Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.legend_handler import HandlerLine2D as HL
import matplotlib.gridspec as gridspec
from obspy.geodetics.base import kilometer2degrees
import seaborn
seaborn.set_style("darkgrid", {"axes.facecolor": ".85"})
#----------------------------------------------------------

def getamp(rfdata, tarray, t):
    """
    Get the amplitude of the receiver function at a specific time

    :param rfdata: receiver function data
    :type rfdata: numpy.ndarray
    :param tarray: time array
    :type tarray: numpy.ndarray
    :param t: time to get the amplitude
    :type t: float

    Returns: 
    Amplitude of the receiver function at time t
    """
    amp = rfdata[(np.abs(tarray - t).argmin())]
    return amp
#----------------------------------------------------------
def hkSeq(rfstreamSed, rfstreamMoho, preonset, HSed=np.linspace(0,10,201),
          KSed=np.linspace(1.65,2.25,201), VpSed=3.0, w1Sed=0.6, w2Sed=0.3, 
          w3Sed=0.1, HMoho=np.linspace(20,60,201), KMoho=np.linspace(1.65,1.95,121), 
          VpMoho=6.9, w1Moho=0.6, w2Moho=0.3, w3Moho=0.1,
          g = [75.,10., 15., 2.5], stack=None):
    """
    Sequential H-K stacking method for receiver function analysis 
    after Yeck et al., 2013
    
    :param rfstreamSed: Stream object of receiver function for sediment layer
                        (high frequency)
    :type rfstreamSed: obspy.core.stream.Stream
    :param rfstreamMoho: Stream object of receiver function for Moho layer
                        (lower frequency)
    :type rfstreamMoho: obspy.core.stream.Stream
    :param preonset: time in seconds before the P-arrival
    :type preonset: integer
    :param HSed: Depth range for sediment layer
    :type HSed: numpy.ndarray
    :param KSed: Vp/Vs range for sediment layer
    :type KSed: numpy.ndarray
    :param VpSed: Vp for sediment layer
    :type VpSed: float
    :param w1Sed: Weight for Ps arrival
    :param w2Sed: Weight for PpPs arrival
    :param w3Sed: Weight for PsPs arrival
    :type w1Sed,w2Sed,w3Sed: float
    :param HMoho: Depth range for Moho layer
    :type HMoho: numpy.ndarray
    :param KMoho: Vp/Vs range for Moho layer
    :type KMoho: numpy.ndarray
    :param VpMoho: Vp for Moho layer
    :type VpMoho: float
    :param w1Moho: Weight for Ps arrival
    :param w2Moho: Weight for PpPs arrival
    :param w3Moho: Weight for PpSs+PsPs arrival
    :type w1Moho,w2Moho,w3Moho: float
    :param g: Gain for plotting
    :type g: list
    :param rmneg: Remove negative values in the stacking results
    :type rmneg: bool
    :param stack: Stack the receiver functions traces into one before 
                  hk stacking
    :type stack: bool

    Returns:
    SequentialHKResult: Dictionary of results from the Sequential 
    H-K stacking method

    Example
    -------

    Initialize the hkSeq module:
    >>> from rfsed.hkSeqYeck import hkSeq
    Define input data and all the necessary parameters. The input data should 
    be in the form of RFStream object from rf package
    >>> import numpy as np
    >>> from rf.rfstream import read_rf
    >>> rfstreamSed = read_rf('path/to/SedimentRF')
    >>> rfstreamMoho = read_rf('path/to/MohoRF')
    >>> preonset = 10
    >>> HSed = np.linspace(0,10,201)
    >>> KSed = np.linspace(1.65,2.25,201)
    >>> VpSed = 3.0
    >>> w1Sed, w2Sed, w3sed = [0.6, 0.3, 0.1]
    >>> HMoho = np.linspace(20,60,201)
    >>> KMoho = np.linspace(1.65,1.95,121)
    >>> VpMoho = 6.9
    >>> w1Moho, w2Moho, w3Moho = [0.6, 0.3, 0.1]
    >>> g = [75.,10., 15., 2.5]
    >>> stack = False
    Call the hkSeq function
    >>> hkSeq(rfstreamSed, rfstreamMoho, preonset, HSed, KSed, VpSed, w1Sed,
            w2Sed, w3Sed, HMoho, KMoho, VpMoho, w1Moho, w2Moho, w3Moho, g, stack)

    """
    if g is None: g = [75.,10., 15., 2.5]
    if VpSed is None: VpSed = 3.0
    if w1Sed is None: w1Sed = 0.6
    if w2Sed is None: w2Sed = 0.3
    if w3Sed is None: w3Sed = 0.1  
    if stack is None: stack = False 
    StackedRFSed=rfstreamSed.stack()
    RFStackedSed=StackedRFSed[0].data
    StackedRFMoho=rfstreamMoho.stack()
    RFStackedMoho=StackedRFMoho[0].data
    staname = rfstreamSed[0].stats.station
    comp =  rfstreamSed[0].stats.channel
    rfdata = rfstreamSed[0].data
    slowness=[]
    for tr in rfstreamMoho:
        slow=tr.stats['slowness']
        slowness.append(slow)
    AvgSlow = np.mean(slowness)
    delta = rfstreamSed[0].stats.delta
    l = len(rfdata)
    t = np.arange(0, l)
    t = (delta *  t) - preonset
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
            for tr in rfstreamSed:
                trdata = tr.data
                b = tr.stats.starttime - tr.stats.onset
                delta= tr.stats.delta
                rp = kilometer2degrees(1) * tr.stats['slowness']
                if stack == True:
                    trdata = StackedRFSed[0].data
                    rp = kilometer2degrees(1) * AvgSlow
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
    print("Best Sediment depth: ", H1, "Best Sediment Vp/Vs:", K1, 
          "Max stack: ", SSed)
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
            for tr in rfstreamMoho:
                trdata = tr.data
                b = rfstreamMoho[0].stats.starttime - rfstreamMoho[0].stats.onset
                delta= tr.stats.delta
                rp = kilometer2degrees(1) * tr.stats['slowness']
                if stack == True:
                    trdata = StackedRFMoho[0].data
                    rp=kilometer2degrees(1) * AvgSlow 
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
    print("Best Moho depth: ", H2, "Best Moho Vp/Vs:", K2, 
          "Max stack: ", SSed)   
    SequentialHKResult = {'RFSed': rfstreamSed, 'RFSedStack': RFStackedSed,
                          'RFMoho': rfstreamMoho, 'RFMohoStack': RFStackedMoho,
                          'time':t, 'bestmodelSed':bmodSed, 
                          'stackvalueSed':stkSed,'bestmodelMoho':bmodMoho,
                          'stackvalueMoho':stkMoho, 'stack':stack, 'comp':comp, 
                          'AvgSlow':AvgSlow, 'VpSed':VpSed, 'VpMoho':VpMoho, 
                          'HSed':HSed, 'KSed':KSed, 'HMoho':HMoho, 
                          'KMoho':KMoho, 'staname':staname}
    return SequentialHKResult  
#----------------------------------------------------------
def plotSeqhk(SequentialHKResult, savepath, g = [75.,10., 15., 2.5], 
              rmneg = True, format = 'jpg'): 
    """
    Plot the results from the Sequential H-K stacking method

    :param SequentialHKResult: Dictionary of results from the Sequential 
                              H-K stacking method (function: hkSequential)
    :type SequentialHKResult: dict
    :param g: Gain for plotting
    :type g: list
    :param rmneg: Remove negative values in the stacking results
    :type rmneg: bool
    :param format: Format of the figure
    :type format: str
    :param savepath: Path to save the figure
    :type savepath: str

    Returns:
    Plot of the results from the Sequential H-K stacking method

    Example
    -------

    Initialize the plotSeqhk module:
    >>> from rfsed.hkSeqYeck import plotSeqhk
    Define input data (which is the result from the hkSeq function) 
    and other plotting parameters. 
    >>> SequentialHKResult = hkSeq(rfstreamSed, rfstreamMoho, preonset, HSed,
                                KSed, VpSed, w1Sed, w2Sed, w3Sed, HMoho, KMoho,
                                VpMoho, w1Moho, w2Moho, w3Moho, g, stack)
    >>> savepath = 'path/to/save/figure'
    >>> g = [75.,10., 15., 2.5]
    >>> rmneg = True
    >>> format = 'pdf'
    Call the plotSeqhk function
    >>> plotSeqhk(SequentialHKResult, savepath, g, rmneg, format)
    
    """
    if rmneg is None: rmneg = True
    if format is None: format = 'pdf'
    if g is None: g = [75.,10., 15., 2.5]
    staname=SequentialHKResult['staname']
    comp=SequentialHKResult['comp']
    stack=SequentialHKResult['stack']
    VpSed=SequentialHKResult['VpSed']
    HSed=SequentialHKResult['HSed']
    KSed=SequentialHKResult['KSed']
    bmodSed = SequentialHKResult['bestmodelSed']
    H1 = bmodSed[0]
    K1 = bmodSed[1]
    stkSed=SequentialHKResult['stackvalueSed']
    if SequentialHKResult['stack'] == True:
        rfstream=SequentialHKResult['RFSedStack']
        rp = kilometer2degrees(1) * SequentialHKResult['AvgSlow']
        t=SequentialHKResult['time']
        term1= ((K1/VpSed)**2 - rp**2)**0.5
        term2= ((1/VpSed)**2 - rp**2)**0.5
        tpsSed = H1 * (term1 - term2)
        tpppsSed = H1 * (term1 + term2)
        tpspsSed = H1 * 2 * (term1)
    else:
        rfstreamSed=SequentialHKResult['RFSed']
        for tr in rfstreamSed:
            rp = kilometer2degrees(1) * tr.stats['slowness']
            term1= ((K1/VpSed)**2 - rp**2)**0.5
            term2= ((1/VpSed)**2 - rp**2)**0.5
            tpsSed = H1 * (term1 - term2)
            tpppsSed = H1 * (term1 + term2)
            tpspsSed = H1 * 2 * (term1)
            tr.stats['tpsSed'] = tpsSed
            tr.stats['tpppsSed'] = tpppsSed
            tr.stats['tpspsSed'] = tpspsSed
    
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
    label='%s Best Model %.2f km %.2f Vp/Vs'%(staname,bmodSed[0],
                                              bmodSed[1]+0.001))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, 
               mode="expand", borderaxespad=0, 
               handler_map={p1:HL(numpoints=1)})
    ax1.set_ylabel('Vp/Vs')
    ax1.set_xlabel('Depth km')
    ax1.set_xlim(min(HSed), max(HSed))
    ax1.set_ylim(min(KSed), max(KSed))
    #----------------------------------------------------------
    # Plot the receiver function witht the estimated times
    #----------------------------------------------------------
    ax2 = plt.subplot(gs[0, 1]) 
    
    if stack == True:
        baz=180
        rfdata = ((SequentialHKResult['RFSedStack']) * g[0]) + baz
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
        ax2.fill_between(t, baz, rfdata, where=rfdata > baz, facecolor='red',
                         alpha = 0.25)
        ax2.fill_between(t, baz, rfdata, where=rfdata <= baz, facecolor='blue',
                         alpha = 0.25)
        mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
        ampgain = g[1]
        ax2.plot([tpsSed, tpsSed], [baz+ampgain, baz-ampgain], 'g-', lw=1.5,
                 label="Moho")
        ax2.plot([tpppsSed, tpppsSed], [baz+ampgain, baz-ampgain], 'g-',
                 lw=1.5)
        ax2.plot([tpspsSed, tpspsSed], [baz+ampgain, baz-ampgain], 'g-',
                 lw=1.5)
        filename = "%s/%s_%s_%s.%s"%(savepath, staname, comp, "Sediment",
                                     format)
        plt.savefig(filename , format=format, transparent=False, dpi=250,
                    bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
        plt.close('all')
    else:
        rfstreamSed = SequentialHKResult['RFSed']
        gcarcarr = np.zeros(len(rfstreamSed))
        for rf in rfstreamSed:        
            delta = rf.stats.delta
            baz = rf.stats.back_azimuth
            b = rf.stats.starttime - rf.stats.onset
            tpsSed = rf.stats.tpsSed
            tpppsSed = rf.stats.tpppsSed
            tpspsSed = rf.stats.tpspsSed
            rfdata = ((rf.data) * g[0]) + baz
            t = (np.arange(0, len(rfdata), 1) * delta) + b
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
            ax2.plot([tpsSed, tpsSed], [baz+ampgain, baz-ampgain], 'g-', 
                     lw=1.5, label="Moho")
            ax2.plot([tpppsSed, tpppsSed], [baz+ampgain, baz-ampgain], 'g-', 
                     lw=1.5)
            ax2.plot([tpspsSed, tpspsSed], [baz+ampgain, baz-ampgain], 'g-', 
                     lw=1.5)
        #----------------------------------------------------------
        # Highlight specifiec time range on RF plot
        #----------------------------------------------------------
        ax2.axvspan(tpsSed-0.1, tpsSed+0.1, alpha=0.5, color='lightgreen')
        ax2.axvspan(tpppsSed-0.1, tpppsSed+0.1, alpha=0.5, color='lightgreen')
        ax2.axvspan(tpspsSed-0.1, tpspsSed+0.1, alpha=0.5, color='lightgreen')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax2.set_xlabel("Time (sec)")
        ax2.set_ylabel("Back-Azimuth (deg)")    
        #----------------------------------------------------------
        # Plot the receiver function witht the estimated times
        #----------------------------------------------------------
        ax3 = plt.subplot(gs[0, 2])
        gcarcarr = np.zeros(len(rfstreamSed))   
        for rf, j in zip(rfstreamSed, np.arange(len(rfstreamSed))):        
            delta = rf.stats.delta
            gcarc = rf.stats.distance
            gcarcarr[j] = gcarc
            b = rf.stats.starttime - rf.stats.onset
            tpsSed = rf.stats.tpsSed
            tpppsSed = rf.stats.tpppsSed
            tpspsSed = rf.stats.tpspsSed
            rfdata = ((rf.data) * g[2]) + gcarc
            t = (np.arange(0, len(rfdata), 1) * delta) + b    
            major_ticks_x = np.arange(-10, 41, 1)                                              
            minor_ticks_x = np.arange(-10, 41, 1) 
            major_ticks_y = np.arange(30, 110, 10)                                              
            minor_ticks_y = np.arange(30, 111, 5)
            ax3.set_xticks(major_ticks_x)                                                       
            ax3.set_xticks(minor_ticks_x, minor=True)                                           
            ax3.set_yticks(major_ticks_y)                                                       
            ax3.set_yticks(minor_ticks_y, minor=True)                                          
            #----------------------------------------------------------
            # Plot grid on the figure
            #----------------------------------------------------------                                                              
            ax3.set_xlim(-1, 10)
            ax3.set_ylim(25,110)
            ax3.plot(t, rfdata, "k-", lw=0.5)
            ax3.fill_between(t, gcarc, rfdata, where=rfdata > gcarc, 
                             facecolor='red', alpha = 0.25)
            ax3.fill_between(t, gcarc, rfdata, where=rfdata <= gcarc, 
                             facecolor='blue', alpha = 0.25)
            mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
            #----------------------------------------------------------
            # Plot H-K sequential formula times
            #----------------------------------------------------------
            ampgain = g[3]
            ax3.plot([tpsSed, tpsSed], [gcarc+ampgain, gcarc-ampgain], 'g-', 
                     lw=1.5, label="Moho")
            ax3.plot([tpppsSed, tpppsSed], [gcarc+ampgain, gcarc-ampgain], 'g-', 
                     lw=1.5)
            ax3.plot([tpspsSed, tpspsSed], [gcarc+ampgain, gcarc-ampgain], 'g-', 
                     lw=1.5)
        #----------------------------------------------------------
        # Highlight specifiec time range on RF plot
        #----------------------------------------------------------
        ax3.axvspan(tpsSed-0.1, tpsSed+0.1, alpha=0.5, color='lightgreen')
        ax3.axvspan(tpppsSed-0.1, tpppsSed+0.1, alpha=0.5, color='lightgreen')
        ax3.axvspan(tpspsSed-0.1, tpspsSed+0.1, alpha=0.5, color='lightgreen')    

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        ax3.set_xlabel("Time (sec)")
        ax3.set_ylabel("Distance (deg)")
        
        filename = "%s/%s_%s_%s.%s"%(savepath, staname, comp, "Sediment", 
                                     format)
        plt.savefig(filename , format=format, transparent=False, dpi=250,
                    bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
        plt.close('all')
    #----------------------------------------------------------
    # Plot Moho Stack
    #----------------------------------------------------------
    VpMoho=SequentialHKResult['VpMoho']
    HMoho=SequentialHKResult['HMoho']
    KMoho=SequentialHKResult['KMoho']
    bmodMoho = SequentialHKResult['bestmodelMoho']  
    H2 = bmodMoho[0]
    K2 = bmodMoho[1]
    stkMoho=SequentialHKResult['stackvalueMoho']

    if SequentialHKResult['stack'] == True:
        rfstream=SequentialHKResult['RFMohoStack']
        rp = kilometer2degrees(1) * SequentialHKResult['AvgSlow']
        t=SequentialHKResult['time']
        term1 = ((K1/VpSed)**2-(rp)**2)**0.5    
        term2 = ((1/VpSed)**2-(rp)**2)**0.5
        term3 = ((K2/VpMoho)**2-(rp)**2)**0.5    
        term4 = ((1/VpMoho)**2-(rp)**2)**0.5
        #--------------------------------------
        tpsMoho = (H1 * (term1 - term2)) + (H2 * (term3 - term4))
        tpppsMoho = (H1 * (term1 + term2)) + (H2 * (term3 + term4))
        tpspsMoho = (2 * H1 * term1) + (2 * H2 * term3)
    else:
        rfstreamMoho=SequentialHKResult['RFMoho']
        for tr in rfstreamMoho:
            rp = kilometer2degrees(1) * tr.stats['slowness']
            #--------------------------------------
            term1 = ((K1/VpSed)**2-(rp)**2)**0.5    
            term2 = ((1/VpSed)**2-(rp)**2)**0.5
            term3 = ((K2/VpMoho)**2-(rp)**2)**0.5    
            term4 = ((1/VpMoho)**2-(rp)**2)**0.5
            #--------------------------------------
            tpsMoho = (H1 * (term1 - term2)) + (H2 * (term3 - term4))
            tpppsMoho = (H1 * (term1 + term2)) + (H2 * (term3 + term4))
            tpspsMoho = (2 * H1 * term1) + (2 * H2 * term3)
            #--------------------------------------      
            tr.stats['tpsMoho'] = tpsMoho
            tr.stats['tpppsMoho'] = tpppsMoho
            tr.stats['tpspsMoho'] = tpspsMoho
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
    label='%s Best Model %.2f km %.2f Vp/Vs'%(staname, bmodMoho[0], 
                                              bmodMoho[1]+0.001))
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
    if stack == True:
        baz=180
        rfdata = ((SequentialHKResult['RFMohoStack']) * g[0]) + baz
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
        ax2.plot([tpsMoho, tpsMoho], [baz+ampgain, baz-ampgain], 'g-', 
                 lw=1.5, label="Moho")
        ax2.plot([tpppsMoho, tpppsMoho], [baz+ampgain, baz-ampgain], 'g-', 
                 lw=1.5)
        ax2.plot([tpspsMoho, tpspsMoho], [baz+ampgain, baz-ampgain], 'g-', 
                 lw=1.5)
        filename = "%s/%s_%s_%s.%s"%(savepath, staname, comp, 'Moho', format)
        plt.savefig(filename , format=format, transparent=False,\
        dpi=250, bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
        plt.close('all')
    else:
        rfstreamMoho= SequentialHKResult['RFMoho']
        gcarcarr = np.zeros(len(rfstreamMoho))
        for rf in rfstreamMoho:        
            delta = rf.stats.delta
            baz = rf.stats.back_azimuth
            b = rf.stats.starttime - rf.stats.onset
            tpsMoho = rf.stats.tpsMoho
            tpppsMoho = rf.stats.tpppsMoho
            tpspsMoho = rf.stats.tpspsMoho
            rfdata = ((rf.data) * g[0]) + baz
            t = (np.arange(0, len(rfdata), 1) * delta) + b
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
            ax2.plot([tpsMoho, tpsMoho], [baz+ampgain, baz-ampgain], 'g-', 
                     lw=1.5, label="Moho")
            ax2.plot([tpppsMoho, tpppsMoho], [baz+ampgain, baz-ampgain], 'g-', 
                     lw=1.5)
            ax2.plot([tpspsMoho, tpspsMoho], [baz+ampgain, baz-ampgain], 'g-', 
                     lw=1.5)
        #----------------------------------------------------------
        # Highlight specifiec time range on RF plot
        #----------------------------------------------------------
        ax2.axvspan(tpsMoho-0.5, tpsMoho+0.5, alpha=0.5, color='lightgreen')
        ax2.axvspan(tpppsMoho-0.5, tpppsMoho+0.5, alpha=0.5, color='lightgreen')
        ax2.axvspan(tpspsMoho-0.5, tpspsMoho+0.5, alpha=0.5, color='lightgreen')
        ax2.set_xlabel("Time (sec)")
        ax2.set_ylabel("Back-Azimuth (deg)")    
        #----------------------------------------------------------
        # Plot the receiver function witht the estimated times
        #----------------------------------------------------------
        ax3 = plt.subplot(gs[0, 2])
        gcarcarr = np.zeros(len(rfstreamMoho))   
        for rf, j in zip(rfstreamMoho, np.arange(len(rfstreamMoho))):        
            delta = rf.stats.delta
            gcarc = rf.stats.distance
            gcarcarr[j] = gcarc
            b = rf.stats.starttime - rf.stats.onset
            tpsMoho = rf.stats.tpsMoho
            tpppsMoho = rf.stats.tpppsMoho
            tpspsMoho = rf.stats.tpspsMoho
            rfdata = ((rf.data) * g[2]) + gcarc
            t = (np.arange(0, len(rfdata), 1) * delta) + b    
            major_ticks_x = np.arange(-10, 41, 5)                                             
            minor_ticks_x = np.arange(-10, 41, 1) 
            major_ticks_y = np.arange(30, 110, 10)                                              
            minor_ticks_y = np.arange(30, 111, 5)
            ax3.set_xticks(major_ticks_x)                                                       
            ax3.set_xticks(minor_ticks_x, minor=True)                                           
            ax3.set_yticks(major_ticks_y)                                                       
            ax3.set_yticks(minor_ticks_y, minor=True)                                          
            #----------------------------------------------------------
            # Plot grid on the figure
            #----------------------------------------------------------                                                               
            ax3.set_xlim(-5, 30)
            ax3.set_ylim(25,110)
            ax3.plot(t, rfdata, "k-", lw=0.5)
            ax3.fill_between(t, gcarc, rfdata, where=rfdata > gcarc, 
                                facecolor='red', alpha = 0.25)
            ax3.fill_between(t, gcarc, rfdata, where=rfdata <= gcarc, 
                                facecolor='blue', alpha = 0.25)
            mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
            #----------------------------------------------------------
            # Plot H-K sequential formula times
            #----------------------------------------------------------
            ampgain = g[3]
            ax3.plot([tpsMoho, tpsMoho], [gcarc+ampgain, gcarc-ampgain], 'g-', 
                     lw=1.5, label="Moho")
            ax3.plot([tpppsMoho, tpppsMoho], [gcarc+ampgain, gcarc-ampgain], 'g-',
                     lw=1.5)
            ax3.plot([tpspsMoho, tpspsMoho], [gcarc+ampgain, gcarc-ampgain], 'g-',
                     lw=1.5)
        #----------------------------------------------------------
        # Highlight specifiec time range on RF plot
        #----------------------------------------------------------
        ax3.axvspan(tpsMoho-0.5, tpsMoho+0.5, alpha=0.5, color='lightgreen')
        ax3.axvspan(tpppsMoho-0.5, tpppsMoho+0.5, alpha=0.5, color='lightgreen')
        ax3.axvspan(tpspsMoho-0.5, tpspsMoho+0.5, alpha=0.5, color='lightgreen')        
        ax3.set_xlabel("Time (sec)")
        ax3.set_ylabel("Distance (deg)")
        filename = "%s/%s_%s_%s.%s"%(savepath, staname, comp, 'Moho', format)
        plt.savefig(filename , format=format, transparent=False, dpi=250, 
                    bbox_inches = 'tight', pad_inches=0.1)
        plt.show()