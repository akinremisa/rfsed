"""
.. module:: ExtractEq
        :synopsis: Extract the data for each earthquake from the local data files
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""

from obspy import read, read_inventory, Stream, read_events, UTCDateTime as UTC
import numpy as np
from obspy.core.event.catalog import read_events
from glob import glob
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from obspy.taup import TauPyModel
#------------------------------------------------------------------------------------
def ExtractEq(datapath, filename, catalog, stalat, stalon, Request_window):
        """
        Extract the data for each earthquake from the local data files
        Writes the extracted data to a new file

        :param datapath: Path to the local data files directory
        :type datapath:str
        :param filename: Path to the new file
        :type filename:str
        :param catalog: Path to the catalog file
        :type catalog:str
        :param stalat: Latitude of the station
        :type stalat:float
        :param stalon: Longitude of the station
        :type stalon:float
        :param Request_window: Time window relative to first P arrival
        :type Request_window:list[start, end]
        :return: Saves extracted earthquake data into a new file
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
                onsettime = i.origins[0].time
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
                rayparam = kilometer2degrees(1) * traveltime[0].ray_param_sec_degree
                incidangle = traveltime[0].incident_angle
                eqstart = (onsettime + P) - abs(Request_window[0])
                eqend = (onsettime + P) + abs(Request_window[1])
                time_window.append([eqstart, eqend])
        #------------------------------------------------------------------------------------
        # Read part of the data with chosen Earthquakes and save to a new file
        Datast = Stream()
        for i in datafiles:
                st = read(i)
                filestart=st[0].stats.starttime
                fileend = st[0].stats.endtime
                for t in time_window:
                        if UTC(t[0]) >= filestart and UTC(t[1]) <= fileend:
                                Eqdata=read(i, starttime=UTC(t[0]), endtime = UTC(t[1]))
                                Datast += Eqdata
        Datast.write(filename)  
#--------------------------------------------------------------------------------