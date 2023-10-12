"""
.. module:: ExtractEqMultiproc
        :synopsis: Extract the data for each earthquake from the local data files in parallel
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""
#%% 
from obspy import read, read_events, Stream, UTCDateTime as UTC
import numpy as np
from obspy.core.event.catalog import read_events
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from obspy.taup import TauPyModel
import multiprocessing
import itertools
def Add_sta_time_2cat(catalogfile, stalat, stalon, Request_window):
    """
    Add station latitude, longitude and request time window to each event in the catalog

    :param catalogfile: Path to the catalog file
    :type catalogfile: str
    :param stalat: Station latitude
    :type stalat: float
    :param stalon: Station longitude
    :type stalon: float
    :param Request_window: Request time window relative to first P arrival (time before P arrival, time after P arrival)
    :type Request_window: list [start, end]
    :return: Catalog with station latitude, longitude and request time window
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
    Get the earthquake start and end time for the requested time window relative to the 
    first P arrival for each event in the catalog. The theoretical P wave arrival time is calculated using the IASP91 reference model

    :param catalog: Catalog with station latitude, longitude and request time window (output of Add_cat_sta_time)
    :type catalog: obspy.core.event.catalog.Catalog
    :return: Dictionary of Earthquake start and end time for the requested time window relative to the first P arrival
    """
    stalat = catalog.stalat
    stalon = catalog.stalon
    Request_window = [catalog.Request_start, catalog.Request_end]       
    onsettime = catalog.origins[0].time
    eqlat = catalog.origins[0].latitude
    eqlon = catalog.origins[0].longitude
    eqmag = catalog.magnitudes[0].mag
    eqdepth = (catalog.origins[0].depth)/1000
    distm,baz,az  = gps2dist_azimuth(stalat, stalon, eqlat, eqlon)
    distkm = distm/1000
    distdeg = locations2degrees(stalat, stalon, eqlat, eqlon)    
    model = TauPyModel(model="iasp91")   
    arrivals = model.get_travel_times(source_depth_in_km=eqdepth, distance_in_degree=distdeg)
    P = arrivals[0].time   
    rayparam = kilometer2degrees(1) * arrivals[0].ray_param_sec_degree
    incidangle = arrivals[0].incident_angle
    eqstart = (onsettime + P) - abs(Request_window[0])
    eqend = (onsettime + P) + abs(Request_window[1])
    timewindow=[{'eqstart':eqstart, 'eqend':eqend}]
    return timewindow
#----------------------------------------------
# Run Parallel processing
def Get_Eqtimes_multiproc(catalog, nproc):
    """
    Run the Eq_times function in parallel for each event in the catalog

    :param catalog: Catalog with station latitude, longitude and request time window (output of Add_cat_sta_time)
    :type catalog: obspy.core.event.catalog.Catalog
    :param nproc: Number of processors to use
    :type nproc: int
    :return: List of dictionaries of Earthquake start and end time for the requested time window relative to the first P arrival
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
    Extract the earthquake data from the local data files for the requested time window relative to the first P arrival

    :param inputparams: List of input parameters (datafiles, timewindow)
    :type inputparams: list
    :return: Stream of earthquake data for the requested time window relative to the first P arrival
    """
    datafiles, timewindow = inputparams
    Datast = Stream()
    Allstream = Stream()
    datafiles = datafiles[0]
    for i in timewindow:
        eqstart=(i[0]['eqstart'])
        eqend=(i[0]['eqend'])
        Eqdata=read(datafiles, starttime=eqstart, endtime = eqend)
        for i in Eqdata:
            Datast.append(i)
    return Datast
# ---------------------------------------------- 
def ExtractEq_Multiproc(datafiles, timewindow, nproc, filename):
    """
    Run the Get_EqStream_multiproc function in parallel for each data file, extracting the earthquake data 
    for the requested time window relative to the first P arrival

    :param datafiles: List of data files
    :type datafiles: list
    :param timewindow: Dictionary of Earthquake start and end time for the requested time window relative to the first P arrival
    :type timewindow: dict
    :param nproc: Number of processors to use
    :type nproc: int
    :param filename: Path to save the extracted earthquake data
    :type filename: str
    :return: Saves the extracted earthquake data to the specified path
    """
    nproc = nproc
    with multiprocessing.Pool(nproc) as pool:
        inputparams=((inputparams, timewindow) for inputparams in itertools.product(datafiles))
        EqStreams=pool.map(Get_EqStream_multiproc, inputparams)   
    pool.close()
    pool.join()
    #----------------------------------------------
    AllStream = Stream()
    for i in range(0,len(EqStreams)):
            AllStream += EqStreams[i]
    AllStream.write(filename)
#----------------------------------------------------------------------------------------------------
