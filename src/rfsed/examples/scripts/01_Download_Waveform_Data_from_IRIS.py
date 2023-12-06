#%% Download data 
import os 
import scipy
import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rfsed.rf import read_rf, RFStream
from rfsed.rf  import get_profile_boxes, iter_event_data, IterMultipleComponents
from rfsed.rf.imaging import plot_profile_map
from rfsed.rf.profile import profile
from tqdm import tqdm
from os.path import exists
from os import mkdir
from rfsed.util import save_Eq_data, save_IRIS_waveform
#------------------------------------------#
savedir=save_IRIS_waveform()
invfile = savedir + '00_eq_stations.xml'
catfile = savedir + '00_eq_events.xml'
datafile = savedir + '00_eq_data.h5'
stname = input("Input Station Name:")
staname= str(stname)
Network = str(input("Input Network Name:"))
StartDate = str(input("Input Start Date (YYYY-MM-DD):"))
EndDate = str(input("Input End Date (YYYY-MM-DD):"))
StationClient=str(input("Input Station Client e.g. IRIS, ORFEUS, GFZ:"))
EventClient=str(input("Input Event Client e.g. IRIS, ORFEUS, GFZ:"))
DataClient=str(input("Input Data Client e.g. IRIS, ORFEUS, GFZ:"))
#------------------------------------------#
# Input Examples
# staname = 'HGN'
# Network = 'NL'
# StartDate = '2015-01-01'
# EndDate = '2015-08-01'
# StationClient='ORFEUS'
# EventClient='GFZ'
# DataClient='ORFEUS'
#------------------------------------------#
# Download station inventory
minlat=50.3; maxlat=53.8  #South-North
minlong=3.0; maxlong=7.7 #West-East
if not os.path.exists(invfile):
    client = Client(StationClient)
    inventory = client.get_stations(network=Network, channel='BH?', level='channel', station=staname,
                                    minlatitude = minlat, maxlatitude = maxlat, 
                                    minlongitude = minlong, maxlongitude = maxlong)
    inventory.write(invfile, 'STATIONXML')
inventory = read_inventory(invfile)
inventory.plot(label=False)
fig = inventory.plot('local')
#------------------------------------------#
# Download event catalog
stacomp= Network + '.' + staname + '..BHZ'
coords = inventory.get_coordinates(stacomp)
lonlat = (coords['longitude'], coords['latitude'])
long=lonlat[0]
lat=lonlat[1]
if not os.path.exists(catfile):
    client = Client(EventClient)
    kwargs = {'starttime': UTC(StartDate), 'endtime': UTC(EndDate),
              'latitude': lonlat[1], 'longitude': lonlat[0],
              'minradius': 30, 'maxradius': 95,
              'minmagnitude': 5.5, 'maxmagnitude': 10}
    catalog = client.get_events(**kwargs)

    catalog.write(catfile, 'QUAKEML')
catalog = read_events(catfile)
fig = catalog.plot(label=None) 
#------------------------------------------#
# Get the data 
if not os.path.exists(datafile):
    client = Client(DataClient)
    stream = RFStream()
    with tqdm() as pbar:
        for s in iter_event_data(catalog, inventory, client.get_waveforms, pbar=pbar):
            stream.extend(s)
    stream.write(datafile, 'H5')  
#%%
eqdata = read_rf(datafile)
for i in eqdata: 
    print(i)
print("Number of Events:", len(eqdata), 'events')
# %%
