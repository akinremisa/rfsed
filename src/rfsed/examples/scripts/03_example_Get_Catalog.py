"""
Extracts the earthquake catalog from IRIS for a given station and time period based on obspy 
"""
# Download the catalog
import os
from os import mkdir
from os.path import exists, dirname
from obspy.clients.fdsn.client import Client
from obspy import UTCDateTime
from obspy.core.event.catalog import read_events
from rfsed.util import save_catalog
savepath=save_catalog()
if not exists(savepath):
    mkdir(savepath)
#Station Information
staname='NE301'
stalat=6.183432
stalon=53.48334
client = Client("IRIS")
starttime = UTCDateTime("2022-02-01")
endtime = UTCDateTime("2022-03-31")
f = "%seq_events.quakeml"%(savepath)
#% Get the catalog and save it
cat = client.get_events(starttime = starttime, endtime = endtime, latitude=stalat, longitude=stalon,\
                        minmagnitude = 6.0, maxmagnitude=10.0, minradius=30, maxradius=90, filename = f)

# Read the file again and plot it
# cat = read_events(f, format="QUAKEML")
# cat.plot()
# %%
