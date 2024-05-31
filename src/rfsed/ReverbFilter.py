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
# module:: ReverbFilter
#      :synopsis: Resonance Filter for receiver function data 
#                 (Removes sediment reverberation effect)
# moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> 
#                Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
# 

import numpy as np
from scipy.signal import detrend
from scipy.fft import fft, ifft
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

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
def ResonanceFilt(rfstream, preonset):
    """
    Resonance Filter for receiver function data 
    (Removes sediment reverberation effect)

    :param rfstream: receiver function data
    :type rfstream: obspy stream
    :param preonset: time in seconds before the P-arrival
    :type preonset: integer

    Returns: 
    A dictionary of stacked receiver function, Filtered receiver function, 
    autocorrelation of the data, resonance filter, time lag (2 way travel time 
    of the sediment reverbration), and the strength of the sediment 
    reverbration (r)

    Example
    -------

    Initialize the ResonanceFilt module:
    >>> from rfsed.ReverbFilter import ResonanceFilt
    Define all the necessary parameters
    The rfstream data is a RFSream object from rf
    >>> from rf.rfstream import read_rf
    >>> rfstream = read_rf('path/to/rfdata')
    >>> preonset = 10
    Call the ResonanceFilt function
    >>> ResonanceFilt(rfstream, preonset)
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
    #---------------------------------------------------------------------
    FltResults={'rf':StackedData, 'filteredrf': filtered_rf, 
                'resonancefilter':resonanceflt, 'time':time, 
                'autoc':ac, 'r':r0, 'tlag':tlag, 'staname':staname, 
                'delta':delta}
    return FltResults
#----------------------------------------------------------

def plotfiltrf(FilteredRF, savepath, format = 'jpg'):
    """
    Plot the filtered receiver function, autocorrelation, and the
    resonance filter

    :param FilteredRF: Dictionary of results from the resonance filtering 
                       method (function: Resonance_Filt)
    :type FilteredRF: dict
    :param savepath: path to save the plots
    :type savepath: str
    :param format: format of the plot (default: jpg)
    :type format: str

    Returns:
    Plots of the filtered receiver function, autocorrelation, and the 
    resonance filter

    Example
    -------

    Initialize the plotfiltrf module:
    >>> from rfsed.ReverbFilter import plotfiltrf
    Define all the necessary parameters
    The FilteredRF data is a dictionary from the ResonanceFilt function
    >>> from rfsed.ReverbFilter import ResonanceFilt
    >>> FilteredRF = ResonanceFilt(rfstream, preonset)
    >>> savepath = 'path/to/saveplots'
    >>> format = 'jpg'
    Call the plotfiltrf function
    >>> plotfiltrf(FilteredRF, savepath, format)
    """
    
    rf=FilteredRF['rf']
    filtered_rf=FilteredRF['filteredrf']
    resonanceflt=FilteredRF['resonancefilter']
    resonanceflt=np.real(resonanceflt)
    time=FilteredRF['time']
    autoc=FilteredRF['autoc']
    staname=FilteredRF['staname']
    # Plot Autocorrelation
    plt.plot(time,autoc)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Autocorrelation of %s'%staname)
    filename = "%s/%s.%s"%(savepath, 'Autocorrelation', format)
    plt.savefig(filename , format=format, transparent=False, dpi=250, 
                bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
    # Plot Resonance Filter
    plt.plot(time, resonanceflt)
    plt.xlabel('Time (s)')
    # plt.ylabel('Amplitude')
    plt.title('Resonance Filter of %s'%staname)
    filename = "%s/%s.%s"%(savepath, "Resonance_Filter", format)
    plt.savefig(filename , format=format, transparent=False, dpi=250,
                bbox_inches = 'tight', pad_inches=0.1)
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
    plt.savefig(filename , format=format, transparent=False, dpi=250,
                bbox_inches = 'tight', pad_inches=0.1)
    plt.show()
    plt.close("all")
#----------------------------------------------------------


