"""
.. module:: ReverbFilter
        :synopsis: Resonance Filter for receiver function data (Removes sediment reverberation effect)
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""

import numpy as np
from scipy.signal import detrend
from scipy.fft import fft, ifft
from scipy.signal import find_peaks
from scipy import signal
from scipy.signal import argrelextrema, argrelmin, argrelmax
import math
import cmath
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
    time=t
    dt=rfstream[0].stats.delta
    StackedRF=rfstream.stack()
    StackedData=StackedRF[0].data 
    #---------------------------------------------
    # Build the frequency vector
    #---------------------------------------------
    n = len(StackedData)
    fmax = 1 / (2.0 * dt)
    df = fmax / (n / 2)
    f = np.hstack((df * np.arange(0, n//2 +1), df * np.arange(-n//2 + 1, 0)))
    nf = n // 2 + 1
    dw = 2.0 * np.pi * df
    w = dw * np.hstack((np.arange(0, n//2+1), np.arange(-n//2 + 1, 0)))
    filtered_rf = np.zeros_like(StackedData)
    D = StackedData
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
    FltResults={'rf':StackedData, 'filteredrf': filtered_rf, 'resonancefilter':resonanceflt, 'time':time, 'autoc':ac, 'r':r0, 'tlag':tlag, 'staname':staname, 'delta':delta}
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


