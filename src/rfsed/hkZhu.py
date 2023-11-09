"""
.. module:: hkZhu
        :synopsis: H-K stacking for receiver functions after Zhu and Kanamori (2000)
        Zhu, L., & Kanamori, H. (2000). Moho depth variation in southern california from
        teleseismic receiver functions. Journal of Geophysical Research: Solid Earth, 105 , 
        2969-2980. doi: 10.1029/1999jb900322
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
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
    :param tarray: time array
    :param t: time to get the amplitude
    :return: amplitude of the receiver function at time t
    """
    amp = rfdata[(np.abs(tarray - t).argmin())]
    return amp
#----------------------------------------------------------

def hk(rfstream, H=np.linspace(20,60,201), K=np.linspace(1.65,1.95,121), Vp=6.9, 
              w1=0.6, w2=0.3, w3=0.1, layer = None, stack=False):
    """
    Calculate the H-K stacking for a receiver function after Zhu and Kanamori (2000)

    :param rfstream: receiver function stream
    :type rfstream: obspy.core.stream.Stream
    :param H: array of Depth values
    :type H: numpy array
    :param K: array of Vp/Vs values
    :type K: numpy array
    :param Vp: P wave velocity of the layer
    :type Vp: float
    :param w1: weight for the Ps arrival
    :param w2: weight for the PpPs arrival
    :param w3: weight for the PsPs+PpSs arrival
    :type w1, w2, w3: float
    :param layer: layer name ('Moho' or 'Sed')
    :type layer: str
    :param stack: stack the receiver functions before the H-K stacking (default False)
    :type stack: bool
    :return: Dictionary of H-K stacking results
    """
    if Vp is None: Vp = 6.5
    if w1 is None: w1 = 0.6
    if w2 is None: w2 = 0.3
    if w3 is None: w3 = 0.1   
    if layer is None: layer = 'Moho'
    if stack is None: stack = False
    StackedRF=rfstream.stack()
    RFStacked=StackedRF[0].data
    slowness=[]
    for tr in rfstream:
        slow=tr.stats['slowness']
        slowness.append(slow)
    AvgSlow = np.mean(slowness)
    rayp=kilometer2degrees(1) * AvgSlow
    #----------------------------------------------------------
    staname = rfstream[0].stats.station
    comp =  rfstream[0].stats.channel
    rfdata = rfstream[0].data
    delta = rfstream[0].stats.delta
    l = len(rfdata)
    t = np.arange(0, l)
    t = (delta *  t) - 10 
    stk = np.zeros((len(K)*len(H),3))
    dbgvalues = np.zeros((len(K)*len(H)*len(rfstream),9))
    z = 0 ; q =0
    for i in range(len(K)):
        Ktemp = K[i]
        for j in range(len(H)):
            Htemp = H[j]
            s = 0.0
            for tr in rfstream:
                trdata = tr.data
                b = tr.stats.starttime - tr.stats.onset
                delta= tr.stats.delta
                rp = kilometer2degrees(1) * tr.stats['slowness']
                if stack == True:
                    trdata = StackedRF[0].data
                    rp = kilometer2degrees(1) * AvgSlow
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
                dbgvalues[q, :] =[rp, Htemp, Ktemp, tps, tppps, tpsps,\
                ampps, ampppps, amppsps]
                q =q + 1
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
    HKResult = {'RF': rfstream, 'RFStack': RFStacked, 'time':t, 'bestmodel':bmod, 'staname':staname, 'Vp':Vp,\
                'stackvalues':stk, 'layer':layer, 'stack':stack, 'comp':comp, 'AvgSlow':AvgSlow, 'H':H, 'K':K}
    return HKResult
#----------------------------------------------------------
def plothk(HKResult, savepath, g = [75.,10., 15., 2.5], rmneg = None, format = 'jpg'): 
    """
    Plot the H-K stacking results

    :param HKResult: Dictionary of H-K stacking results from hk function
    :type HKResult: dict
    :param savepath: path to save the figure
    :type savepath: str
    :param g: list of gain values for the plot
    :type g: list
    :param rmneg: remove the negative values from the stack (default True)
    :type rmneg: bool
    :param format: format of the figure
    :type format: str
    """
    if rmneg is None: rmneg = True
    if format is None: format = 'pdf'
    if g is None: g = [75.,10., 15., 2.5]
    staname=HKResult['staname']
    bmod = HKResult['bestmodel']
    Hbm = bmod[0]
    Kbm = bmod[1]
    Sbm = bmod[2]
    layer = HKResult['layer']
    comp=HKResult['comp']
    Vp=HKResult['Vp']
    if HKResult['stack'] == True:
        rfstream=HKResult['RFStack']
        rp = kilometer2degrees(1) * HKResult['AvgSlow']
        t=HKResult['time']
        term1= ((Kbm/Vp)**2 - rp**2)**0.5
        term2= ((1/Vp)**2 - rp**2)**0.5
        tps = Hbm * (term1 - term2)
        tppps = Hbm * (term1 + term2)
        tpsps = Hbm * 2 * (term1)
    else:
        rfstream=HKResult['RF']
        for tr in rfstream:
            rp = kilometer2degrees(1) * tr.stats['slowness']
            term1= ((Kbm/Vp)**2 - rp**2)**0.5
            term2= ((1/Vp)**2 - rp**2)**0.5
            tps = Hbm * (term1 - term2)
            tppps = Hbm * (term1 + term2)
            tpsps = Hbm * 2 * (term1)    
            tr.stats['tps'] = tps
            tr.stats['tppps'] = tppps
            tr.stats['tpsps'] = tpsps
    #----------------------------------------------------------
    # plot H-K results
    #----------------------------------------------------------
    if HKResult['layer'] == 'Moho':
        mintime=-5
        maxtime=30
        span=0.5
        timestep=5
    elif HKResult['layer'] == 'Sed':
        mintime=-1
        maxtime=10
        span=0.1
        timestep=1
    H=HKResult['H']
    K=HKResult['K']
    stk = HKResult['stackvalues']
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
    if HKResult['stack'] == True:
        rfstream=HKResult['RFStack']
        baz=180
        rfdata =  ((rfstream) * g[0]) + baz
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
        ax2.axvspan(tps-span, tps+span, alpha=0.5, color='lightgreen')
        ax2.axvspan(tppps-span, tppps+span, alpha=0.5, color='lightgreen')
        ax2.axvspan(tpsps-span, tpsps+span, alpha=0.5, color='lightgreen')
        ax2.set_xlabel("Time (sec)")
        ax2.set_ylabel("Back-Azimuth (deg)") 
        filename = "%s/%s_%s_%s.%s"%(savepath, staname, comp, str(layer), format)
        plt.savefig(filename , format=format, transparent=False,\
        dpi=250, bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
        plt.close('all')
    else:
        rfstream=HKResult['RF']
        gcarcarr = np.zeros(len(rfstream))
        for rf in rfstream:        
            delta = rf.stats.delta
            baz = rf.stats.back_azimuth
            b = rf.stats.starttime - rf.stats.onset
            tps = rf.stats.tps
            tppps = rf.stats.tppps
            tpsps = rf.stats.tpsps
            rfdata = ((rf.data) * g[0]) + baz
            t = (np.arange(0, len(rfdata), 1) * delta) + b
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
        # Plot the receiver function by distance with the estimated times 
        #----------------------------------------------------------
        ax3 = plt.subplot(gs[0, 2])
        gcarcarr = np.zeros(len(rfstream))   
        for rf, j in zip(rfstream, np.arange(len(rfstream))):        
            delta = rf.stats.delta
            gcarc = rf.stats.distance
            gcarcarr[j] = gcarc
            b = rf.stats.starttime - rf.stats.onset
            tps = rf.stats.tps
            tppps = rf.stats.tppps
            tpsps = rf.stats.tpsps
            rfdata = ((rf.data) * g[2]) + gcarc
            t = (np.arange(0, len(rfdata), 1) * delta) + b    
            major_ticks_x = np.arange(-10, 41, timestep)                                              
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
            ax3.set_xlim(mintime, maxtime)
            ax3.set_ylim(25,110)
            ax3.plot(t, rfdata, "k-", lw=0.5)
            ax3.fill_between(t, gcarc, rfdata, where=rfdata > gcarc, facecolor='red', alpha = 0.25)
            ax3.fill_between(t, gcarc, rfdata, where=rfdata <= gcarc, facecolor='blue', alpha = 0.25)
            mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
            #----------------------------------------------------------
            # Plot H-K sequential formula times
            #----------------------------------------------------------
            ampgain = g[3]
            ax3.plot([tps, tps], [gcarc+ampgain, gcarc-ampgain], 'g-', lw=1.5)
            ax3.plot([tppps, tppps], [gcarc+ampgain, gcarc-ampgain], 'g-', lw=1.5)
            ax3.plot([tpsps, tpsps], [gcarc+ampgain, gcarc-ampgain], 'g-', lw=1.5)
        #----------------------------------------------------------
        # Highlight specifiec time range on RF plot
        #----------------------------------------------------------
        ax3.axvspan(tps-span, tps+span, alpha=0.5, color='lightgreen')
        ax3.axvspan(tppps-span, tppps+span, alpha=0.5, color='lightgreen')
        ax3.axvspan(tpsps-span, tpsps+span, alpha=0.5, color='lightgreen')    
        ax3.set_xlabel("Time (sec)")
        ax3.set_ylabel("Distance (deg)")
        filename = "%s/%s_%s_%s.%s"%(savepath, staname, comp, str(layer), format)
        plt.savefig(filename , format=format, transparent=False,\
            dpi=250, bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
