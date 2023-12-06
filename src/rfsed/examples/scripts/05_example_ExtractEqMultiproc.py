# %% Test the ExtractEqMultiproc Function
from glob import glob
import os
from obspy import read
from rfsed.ExtractEqMultiproc import Add_sta_time_2cat, Get_Eqtimes_multiproc, ExtractEq_Multiproc
from rfsed.util import catalog_example, raw_data_example, save_Eq_data

staname='NE301'
stalat=6.183432
stalon=53.48334
datapath= raw_data_example()
datafiles = sorted(glob("%s*.dat"%(datapath)))
savedir=save_Eq_data()
savepath=savedir+'/Extracted_Eq_Data_Multiproc/'
if not os.path.exists(savepath):
    os.mkdir(savepath)
filename = savepath + 'eq_data.h5'
catalogfile=catalog_example()
Request_window=[-50, 150] #Time relative to first P arrival
nproc = 40 #Number of processors to use if parallel processing is used
catalog = Add_sta_time_2cat(catalogfile, stalat, stalon, Request_window)
timewindow=Get_Eqtimes_multiproc(catalog, nproc)
eqdata = ExtractEq_Multiproc(datafiles=datafiles, timewindow=timewindow, nproc=nproc, filename=filename)
#%% Check the number of traces in the file
data_extract=read(filename )
print(len(data_extract), 'Traces in the file')
print(data_extract)