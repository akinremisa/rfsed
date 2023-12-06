# %% Run the function to extract the data from local files
import os
from obspy import read
from os.path import exists
from os import mkdir
from rfsed.ExtractEq import ExtractEq
from rfsed.util import catalog_example, raw_data_example, save_Eq_data
staname='NE301'
stalat=6.183432
stalon=53.48334
datapath= raw_data_example()
savedir=save_Eq_data()
savepath=savedir+'/Extracted_Eq_Data/'
if not exists(savepath):
    mkdir(savepath)
filename = savepath + 'eq_data.h5'
catalog = catalog_example()
Request_window=[-50, 150] #Time relative to first P arrival
ExtractData= ExtractEq(datapath, filename, catalog, stalat, stalon, Request_window)
#%% Check the extracted data
ExtractedEQ_Data=read(filename)
print(len(ExtractedEQ_Data), ' traces in the extracted data')
print(ExtractedEQ_Data)
# %%
