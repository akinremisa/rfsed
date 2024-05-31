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

#
# module:: hkYu
#      :synopsis: Modified H-K stacking method of Yu et al. (2015) for 
#      receiver functions filtered with the resonance filter
#      Yu, Y., Song, J., Liu, K. H., & Gao, S. S. (2015). Determining crustal 
#      structure beneath seismic stations overlying a low-velocity sedimentary 
#      layer using receiver functions. Journal of Geophysical Research: Solid 
#      Earth, 120 , 3208-3218. doi:10.1002/2014JB011610
# moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> 
#                Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
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
    
    Returns: 
    Amplitude of the receiver function at time t
    """
    amp = rfdata[(np.abs(tarray - t).argmin())]
    return amp
#----------------------------------------------------------

def hkYu(FltResult, rayp, HSubSed, KSubSed, VpMoho, VpSed, VsSed, w1SubSed, 
         w2SubSed, w3SubSed):
    """
    Modified H-K stacking method of Yu et al. (2015)

    :param FltResult: Dictionary of results from the resonance filtering 
                      method (function: Resonance_Filt)
    :type FltResult: dict
    :param rayp: Ray parameter
    :type rayp: float
    :param HSubSed: Subsediment layer thickness array
    :type HSubSed: numpy array
    :param KSubSed: Subsediment layer Vp/Vs array
    :type KSubSed: numpy array
    :param VpMoho: Moho Vp
    :type VpMoho: float
    :param VpSed: Sediment Vp
    :type VpSed: float
    :param w1SubSed: Weight for the Subsediment Ps arrival at adjusted 
                    arrival time
    :type w1SubSed: float
    :param w2SubSed: Weight for the Subsediment Ppps arrival at adjusted 
                    arrival time
    :type w2SubSed: float
    :param w3SubSed: Weight for the Subsediment Psps+PpSs arrival at adjusted 
                    arrival time
    :type w3SubSed: float

    Returns:
    HKYuResult: Dictionary of results from the mofified H-K stacking method

    Example
    -------

    Initialize the Modified H-K stacking method:
    >>> from rfsed.hkYu import hkYu
    >>> import numpy as np
    Define all the necessary parameters
    The FltResult are results from the resonance filtering method 
    (see ReverbFilter module)
    >>> FltResult = {'filteredrf': filteredrf, 'time':time, 'tlag':tlag}
    >>> rayp = 0.04
    >>> HSubSed = np.linspace(20,60,201)
    >>> KSubSed = np.linspace(1.65,1.95,121)
    >>> VpMoho = 6.9
    >>> VpSed = 2.5
    >>> VsSed = 1.4
    >>> w1SubSed, w2SubSed, w3SubSed = [0.6, 0.3, 0.2]
    Call the hkYu function
    >>> hkYuResult = hkYu(FltResult, rayp, HSubSed, KSubSed, VpMoho, VpSed, 
                            VsSed, w1SubSed, w2SubSed, w3SubSed)

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
    # VsSed=0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
    SedH=tlag*VsSed/2
    print('Sediment thickness based on the layer 2-way travel time is',  
          SedH, 'km' )
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
    plt.tricontourf(stkSubSed[:,0], stkSubSed[:,1], stkSubSed[:,2],60, 
                    cmap='jet')
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    p1,=plt.plot(bmodSubSed[0],bmodSubSed[1], 'k+', mew=5, ms=15, 
                 label='Best Model %s km %s Vp/Vs'%(bmodSubSed[0], 
                                                    bmodSubSed[1]))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, 
               mode="expand", borderaxespad=0,
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
    print("Best Subsediment thickness: ", HSubSedbm, "Best Subsediment Vp/Vs:", 
          KSubSedbm, "Max stack: ", SSubSedbm)
    HKYuResult = {'FltRF': Fltrf, 'time':t, 'tlag':tlag, 'Pdelay':Pdelay, 
                  'rayp':rayp,'bestmodelSubSed':bmodSubSed, 'VpMoho':VpMoho, 
                  'VpSed':VpSed, 'stackvaluesSubSed':stkSubSed, 
                  'HMoho':HSubSed, 'KMoho':KSubSed, 'SubSedThick':HSubSedbm,
                  'SubSedVpVs':KSubSedbm,'SedThick':SedH}
    #----------------------------------------------------------
    return HKYuResult
#----------------------------------------------------------

def plothkYu(hkYuResult, savepath, g = [75.,10., 15., 2.5], rmneg = True, 
             format = 'jpg'): 
    """
    Plot the results of the Modified H-K stacking method of Yu et al. (2015)

    :param hkYuResult: Dictionary of results from the Modified H-K stacking 
                        method (function: hkYu)
    :type hkYuResult: dict
    :param savepath: path to save the plots
    :type savepath: str
    :param g: gain values for the plot (default: [75.,10., 15., 2.5])
    :type g: list
    :param rmneg: remove the stack values that is lower than zero 
                  (default: False)
    :type rmneg: bool
    :param format: format of the plot (default: jpg)
    :type format: str   

    Returns:
    Plot of the results from the Modified H-K stacking method

    Example
    -------

    Initialize the Modified H-K stacking plotting method:
    >>> from rfsed.hkYu import plothkYu
    Define all the necessary parameters
    The hkYuResult are results from the Modified H-K stacking 
    method (see hkYu function)
    >>> hkYuResult = hkYu(FltResult, rayp, HSubSed, KSubSed, HSed, KSed,
                        VpMoho, VpSed, VsSed, w1SubSed, w2SubSed, w3SubSed)
    >>> savepath = 'path/to/saveplots'
    >>> g = [75.,10., 15., 2.5]
    >>> rmneg = True
    >>> format = 'jpg'
    Call the plothkYu function
    >>> plothkYu(hkYuResult, savepath, g, rmneg, format)
    """
    
    if rmneg is None: rmneg = False
    if format is None: format = 'pdf'
    if g is None: g = [75.,10., 15., 2.5]
    staname='FilteredRF_HKYu'
    Pdelay=hkYuResult['Pdelay']
    tlag=hkYuResult['tlag']
    VpMoho=hkYuResult['VpMoho']
    VpSed=hkYuResult['VpSed']
    HMoho=hkYuResult['HMoho']
    KMoho=hkYuResult['KMoho']
    rfdata=hkYuResult['FltRF']
    rp = hkYuResult['rayp']
    t=hkYuResult['time']
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
    plt.tricontourf(stkSubSed[:,0], stkSubSed[:,1], stkSubSed[:,2],50, 
                    cmap ="jet")
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    cb.ax.set_xlabel('S', rotation=0)
    p1,=ax1.plot(bmodSubSed[0],bmodSubSed[1], 'k+', mew=2, ms=10, 
                 label='Best SubSediment Model %.2f km %.2f Vp/Vs'%(bmodSubSed[0],
                                                                    bmodSubSed[1]+0.001))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, 
               mode="expand", borderaxespad=0, 
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
    ax2.plot([tpsSubSed, tpsSubSed], [baz+ampgain, baz-ampgain], 'g-', 
             lw=1.5, label="Moho")
    ax2.plot([tpppsSubSed, tpppsSubSed], [baz+ampgain, baz-ampgain], 'g-', 
             lw=1.5)
    ax2.plot([tpspsSubSed, tpspsSubSed], [baz+ampgain, baz-ampgain], 'g-', 
             lw=1.5)
    filename = "%s/%s_%s.%s"%(savepath, staname, 'SubSediment', format)
    plt.savefig(filename , format=format, transparent=False,\
    dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close('all')
#----------------------------------------------------------
