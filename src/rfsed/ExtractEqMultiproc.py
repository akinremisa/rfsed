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
# module:: ExtractEqMultiproc
#      :synopsis: Extract the data for each earthquake from the local 
#                    data files in parallel
# moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> 
#                   Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
#

from obspy import read, read_events, Stream
import numpy as np
from obspy.core.event.catalog import read_events
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from obspy.taup import TauPyModel
import multiprocessing
import itertools

def Add_sta_time_2cat(catalogfile, stalat, stalon, Request_window):
    """
    Add station latitude, longitude and request time window to each 
    event in the catalog

    :param catalogfile: Path to the catalog file
    :type catalogfile: str
    :param stalat: Station latitude
    :type stalat: float
    :param stalon: Station longitude
    :type stalon: float
    :param Request_window: Request time window relative to first P arrival 
                           (time before P arrival, time after P arrival)
    :type Request_window: list [start, end]


    Returns:
        Catalog of events with station latitude, longitude and request time 
        window added to each event

    Example
    --------

    >>> # Initialize the Add_sta_time_2cat function:
    >>> from rfsed.ExtractEqMultiproc import Add_sta_time_2cat
    >>> # Define the input parameters:
    >>> catalogfile = 'path/to/catalogfile'
    >>> stalat = 52.0
    >>> stalon = 4.0
    >>> Request_window = [-50, 150]
    >>> # Add the station latitude, longitude and request time window to the 
    >>> # catalog using the Add_sta_time_2cat function:
    >>> catalog = Add_sta_time_2cat(catalogfile, stalat, stalon,
                                    Request_window)
    """

    catalog = read_events(catalogfile)
    for i in catalog:
        i.stalat = stalat
        i.stalon = stalon
        i.Request_start = Request_window[0]
        i.Request_end = Request_window[1]
    return catalog
#----------------------------------------------
def Eq_times(catalog):
    """
    Get the earthquake start and end times in obspy UTCDateTime for the 
    requested time window relative to the first P arrival for each event in 
    the catalog. The theoretical P wave arrival time is calculated using the 
    IASP91 reference model

    :param catalog: Catalog with station latitude, longitude and request 
                    time window (output of Add_cat_sta_time)
    :type catalog: obspy.core.event.catalog.Catalog

    Returns: 
        Dictionary of Earthquake start and end time for the requested time window 
        relative to the first P arrival

    Example
    --------

    >>> # Initialize the Eq_times function:
    >>> from rfsed.ExtractEqMultiproc import Eq_times
    >>> # Define the input parameters (output from the Add_cat_sta_time function):
    >>> catalog = Add_cat_sta_time(catalogfile, stalat, stalon, Request_window)
    >>> # Get the earthquake start and end times for the requested time window 
    >>> # relative to the first P arrival using the Eq_times function:
    >>> timewindow = Eq_times(catalog)

    """
    stalat = catalog.stalat
    stalon = catalog.stalon
    Request_window = [catalog.Request_start, catalog.Request_end]       
    eqtime = catalog.origins[0].time
    eqlat = catalog.origins[0].latitude
    eqlon = catalog.origins[0].longitude
    eqmag = catalog.magnitudes[0].mag
    eqdepth = (catalog.origins[0].depth)/1000
    distm,baz,az  = gps2dist_azimuth(stalat, stalon, eqlat, eqlon)
    distkm = distm/1000
    distdeg = locations2degrees(stalat, stalon, eqlat, eqlon)    
    model = TauPyModel(model="iasp91")   
    arrivals = model.get_travel_times(source_depth_in_km=eqdepth, 
                                      distance_in_degree=distdeg)
    P = arrivals[0].time 
    onset=eqtime + P  
    rayparam = kilometer2degrees(1) * arrivals[0].ray_param_sec_degree
    incidangle = arrivals[0].incident_angle
    eqstart = (eqtime + P) - abs(Request_window[0])
    eqend = (eqtime + P) + abs(Request_window[1])
    timewindow=[{'eqstart':eqstart, 'eqend':eqend, 'onset':onset}]
    return timewindow
#----------------------------------------------

# Run Parallel processing
def Get_Eqtimes_multiproc(catalog, nproc):
    """
    Run the Eq_times function in parallel for all events in the catalog

    :param catalog: Catalog with station latitude, longitude and request
                    time window (output of Add_cat_sta_time)
    :type catalog: obspy.core.event.catalog.Catalog
    :param nproc: Number of processors to use
    :type nproc: int

    
    Returns: 
        List of dictionaries of Earthquake start and end times 
        (in obspy UTCDateTime) for the requested time window relative to the first P arrival

    Example
    --------

    >>> # Initialize the Get_Eqtimes_multiproc function:
    >>> from rfsed.ExtractEqMultiproc import Get_Eqtimes_multiproc
    >>> # Define the input parameters catalog and nproc (number of processors to use)
    >>> # he catalog is an output from the Add_cat_sta_time function
    >>> catalog = Add_cat_sta_time(catalogfile, stalat, stalon, Request_window)
    >>> nproc = 4
    >>> # Get the earthquake start and end times for the requested time window 
    >>> # relative to the first P arrival using the Get_Eqtimes_multiproc function:
    >>> timewindow = Get_Eqtimes_multiproc(catalog, nproc)
    """

    nproc = nproc
    pool = multiprocessing.get_context('fork').Pool(nproc)
    timewindow = pool.map(Eq_times, catalog)
    pool.close()
    pool.join()
    return timewindow
#----------------------------------------------

def Get_EqStream_multiproc(inputparams):
    """
    Extract the earthquake waveform data from the local data files for the
    requested time window relative to the first P arrival

    :param inputparams: List of input parameters (datafiles, timewindow)
    :type inputparams: list

    
    Returns:
        Stream of earthquake data for the requested time window relative 
        to the first P arrival

    Example
    --------

    >>> # Initialize the Get_EqStream_multiproc function:
    >>> from rfsed.ExtractEqMultiproc import Get_EqStream_multiproc
    >>> # Define the input parameters datafiles and timewindow
    >>> datafiles = ['path/to/datafiles']
    >>> timewindow = [{'eqstart':eqstart, 'eqend':eqend, 'onset':onset}]
    >>> inputparams = (datafiles, timewindow)
    >>> # Get the earthquake waveform data for the requested time window relative 
    >>> # to the first P arrival using the Get_EqStream_multiproc function:
    >>> Datast = Get_EqStream_multiproc(inputparams)
    """

    datafiles, timewindow = inputparams
    Datast = Stream()
    Allstream = Stream()
    datafiles = datafiles[0]
    for i in timewindow:
        eqstart=(i[0]['eqstart'])
        eqend=(i[0]['eqend'])
        onset=(i[0]['onset'])
        Eqdata=read(datafiles, starttime=eqstart, endtime = eqend)
        for i in Eqdata:
            i.stats.onset=onset
            Datast.append(i)
    return Datast
# ---------------------------------------------- 

def ExtractEq_Multiproc(datafiles, timewindow, nproc, filename):
    """
    Runs the Get_EqStream_multiproc function in parallel for each data file, 
    extracts the earthquake data for the requested time window relative to the 
    first P arrival and save the data to a file

    :param datafiles: List of data files
    :type datafiles: list
    :param timewindow: Dictionary of Earthquake start and end time for the 
    requested time window relative to the first P arrival
    :type timewindow: dict
    :param nproc: Number of processors to use
    :type nproc: int
    :param filename: Path to save the extracted earthquake data
    :type filename: str

    
    Returns: 
        Extracted earthquake data to the specified file path

    Example
    --------

    >>> # Initialize the ExtractEq_Multiproc function:
    >>> from rfsed.ExtractEqMultiproc import ExtractEq_Multiproc
    >>> # Define the input parameters datafiles, timewindow, nproc and filename
    >>> datafiles = ['path/to/datafiles']
    >>> # time window is an output from the Eq_times or Get_Eqtimes_multiproc 
    >>> # functions
    >>> timewindow = [{'eqstart':eqstart, 'eqend':eqend, 'onset':onset}]
    >>> nproc = 4
    >>> filename = 'path/to/newfile'
    >>> # Extract the earthquake waveform data for the requested time window relative
    >>> # to the first P arrival using the ExtractEq_Multiproc function:
    >>> ExtractEq_Multiproc(datafiles, timewindow, nproc, filename)
    """
    
    nproc = nproc
    with multiprocessing.Pool(nproc) as pool:
        inputparams=((inputparams, timewindow) for inputparams in 
                     itertools.product(datafiles))
        EqStreams=pool.map(Get_EqStream_multiproc, inputparams)   
    pool.close()
    pool.join()
    #----------------------------------------------
    AllStream = Stream()
    for i in range(0,len(EqStreams)):
            AllStream += EqStreams[i]
    AllStream.write(filename)
#----------------------------------------------------------------------------------------------------
