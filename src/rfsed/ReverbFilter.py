"""
.. module:: ReverbFilter
        :synopsis: Resonance Filter for receiver function data (Removes sediment reverberation effect)
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
def Resonance_Filt(rfstream, preonset):
    """
    Resonance Filter for receiver function data (Removes sediment reverberation effect)

    :param rfstream: receiver function data
    :type rfstream: obspy stream
    :param preonset: time in seconds before the P-arrival
    :type preonset: integer

    :return: A dictionary of stacked receiver function, Filtered receiver function, 
            autocorrelation of the data, resonance filter, time lag (2 way travel time of the sediment reverbration), 
            and the strength of the sediment reverbration (r)
    """
    delta = rfstream[0].stats.delta
    staname=rfstream[0].stats.station
    comp =  rfstream[0].stats.channel
    l = len(rfstream[0].data)
    t = np.arange(0, l)
    t = (delta *  t) - preonset 
    Dt=rfstream[0].stats.delta
    StackedRF=rfstream.stack()
    StackedData=StackedRF[0].data
    tzero=np.where(t==0)[0]
    tzero=int(tzero)  
    start=tzero
    rf=StackedData[start::] 
    time = t[start::]
    Dt=delta
    N=len(rf)
    fmax=1/(2*Dt)
    df=fmax/(N/2)
    #---------------------------------------------
    # Build the frequency vector
    N2=N/2
    #---------------------------------------------
    column1=np.arange(0,N2+1)    
    column2=np.arange(-N2+1,-1) 
    #---------------------------------------------
    # Join the two vectors
    full_column = np.concatenate((column1,column2))
    #---------------------------------------------
    f=df*full_column
    Nf = N/2+1
    dw=2*np.pi*df
    w=dw*full_column
    #----------------------------------------------------------------------------------------------------------------------
    # Do Auto-correlation
    #----------------------------------------------------------------------------------------------------------------------
    D=rf
    D=D-np.mean(D)
    D=signal.detrend(D, axis=-1, type='linear')
    #------------------------
    autoc=signal.correlate(D,D,mode='full')#,method='auto')
    autoc=autoc/max(autoc)
    autoc=autoc[N-1:2*N-1]
    #----------------------------------------------------------------------------------------------------------------------
    # Find time lag and amplitude of the autocorrelation
    #----------------------------------------------------------
    local_min=argrelmin(autoc)
    index=(local_min[0][0])
    r=(autoc[index])
    r=math.sqrt(r**2)
    tlag= time[index]
    #------------------------------------------------------------------------------
    # Build Filter
    #------------------------------------------------------------------------------
    resonanceflt = []
    for i in range(0, len(w)):
        flt= (1 + (r)*cmath.exp((-1j)*(w[i])*tlag))
        resonanceflt.append(flt)
    #------------------------------------------------------------------------------
    # Filter the Data with a Resonance Filter (Fourier Transform)
    F_trans=fft(rf)
    Filtered=F_trans*resonanceflt
    F_inv=ifft(Filtered)
    #Use Real Part of the Inverse Fourier Transform
    filtered_rf=np.real(F_inv)  #real part of the inverse fourier transform
    #------------------------------------------------------------------------------
    FltResults={'rf':rf, 'filteredrf': filtered_rf, 'resonancefilter':resonanceflt, 'time':time, 'autoc':autoc, 'r':r, 'tlag':tlag, 'staname':staname, 'delta':delta}
    return FltResults
#----------------------------------------------------------
def plotfiltrf(FilteredRF, savepath, format = 'jpg'):
    """
    Plot the filtered receiver function, autocorrelation, and the resonance filter

    :param FilteredRF: Dictionary of results from the resonance filtering method (function: Resonance_Filt)
    :type FilteredRF: dict
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
    staname=FilteredRF['staname']
    # Plot Autocorrelation
    plt.plot(time,autoc)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Autocorrelation of %s'%staname)
    filename = "%s/%s.%s"%(savepath, 'Autocorrelation', format)
    plt.savefig(filename , format=format, transparent=False,\
            dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
    # Plot Resonance Filter
    plt.plot(time, resonanceflt)
    plt.xlabel('Time (s)')
    # plt.ylabel('Amplitude')
    plt.title('Resonance Filter of %s'%staname)
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
    plt.title('RF and Filtered RF of %s'%staname)
    plt.legend()
    filename = "%s/%s.%s"%(savepath, "RF_FilteredRF", format)
    plt.savefig(filename , format=format, transparent=False,\
            dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
#----------------------------------------------------------


