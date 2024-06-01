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
# module:: ExtractEq
#        :synopsis: Extract the earthquake waveform data  from the local data 
#                    files
# moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> 
#                Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
# 

from obspy import read, read_inventory, Stream, read_events, UTCDateTime as UTC
import numpy as np
from obspy.core.event.catalog import read_events
from glob import glob
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from obspy.taup import TauPyModel
#------------------------------------------------------------------------------
def ExtractEq(datapath, filename, catalog, stalat, stalon, Request_window):
    """
    Extract the earthquake waveform data from the local data files and writes 
    the extracted data to a new file

    :param datapath: Path to the local data files directory
    :type datapath: str
    :param filename: Path to the new file
    :type filename: str
    :param catalog: Path to the catalog file
    :type catalog: str
    :param stalat: Latitude of the station
    :type stalat: float
    :param stalon: Longitude of the station
    :type stalon: float
    :param Request_window: Time window relative to first P arrival (in seconds)
    :type Request_window: list[start, end]

    
    Returns:  
        Extracted earthquake data in a new file

    Example
    -------

    >>> # Initialize the ExtractEq module:
    >>> from rfsed.ExtractEq import ExtractEq
    >>> # Define all the necessary parameters
    >>> datapath = 'path/to/datafiles/'
    >>> filename = 'path/to/newfile'
    >>> catalog = 'path/to/catalog'
    >>> stalat = 52.22
    >>> stalon = 6.89
    >>> Request_window = [-50, 150]
    >>> # Call the ExtractEq function
    >>> ExtractEq(datapath, filename, catalog, stalat, stalon, Request_window)

    """
    datafiles = sorted(glob("%s*.dat"%(datapath)))
    len_datafiles = len(datafiles)
    firstfile = read(datafiles[0])
    lastfile = read(datafiles[len_datafiles - 1])
    recordstarttime = firstfile[0].stats.starttime 
    recordendtime = lastfile[0].stats.endtime
    # staname=firstfile[0].stats.station
    # print('Start of all data record = ', recordstarttime)
    # print('End of all data record = ', recordendtime)
    #----------------------------------------------
    # Get the requested UTC time window for each earthquake
    catalog = read_events(catalog)
    len_catalog=len(catalog)
    # print(len_catalog, ' events in Catalog')    
    time_window = []
    for i in catalog:
        eqtime = i.origins[0].time
        eqlat = i.origins[0].latitude
        eqlon = i.origins[0].longitude
        eqmag = i.magnitudes[0].mag
        eqdepth = (i.origins[0].depth)/1000
        distm,baz,az  = gps2dist_azimuth(stalat, stalon, eqlat, eqlon)
        distkm = distm/1000
        distdeg = locations2degrees(stalat, stalon, eqlat, eqlon)  
        model = TauPyModel(model="iasp91")  
        traveltime = model.get_travel_times(source_depth_in_km=eqdepth, distance_in_degree=distdeg)
        arrivals=traveltime
        P = traveltime[0].time  
        onset=eqtime + P 
        rayparam = kilometer2degrees(1) * traveltime[0].ray_param_sec_degree
        incidangle = traveltime[0].incident_angle
        eqstart = (eqtime + P) - abs(Request_window[0])
        eqend = (eqtime + P) + abs(Request_window[1])
        time_window.append([eqstart, eqend, onset])
    #----------------------------------------------------------------------
    # Read part of the data with chosen Earthquakes and save to a new file
    Datast = Stream()
    for i in datafiles:
        st = read(i)
        filestart=st[0].stats.starttime
        fileend = st[0].stats.endtime
        for t in time_window:
                onset = t[2]
                if UTC(t[0]) >= filestart and UTC(t[1]) <= fileend:
                        Eqdata=read(i, starttime=UTC(t[0]), 
                                        endtime = UTC(t[1]))
                        tr=len(Eqdata)
                        Eqdata[0].stats.onset=onset
                        try:
                                Eqdata[1].stats.onset=onset
                                Eqdata[2].stats.onset=onset
                                Eqdata[3].stats.onset=onset
                                Eqdata[4].stats.onset=onset
                                Eqdata[5].stats.onset=onset
                        except:
                                continue
                        Datast += Eqdata
    Datast.write(filename)  
#--------------------------------------------------------------------------------