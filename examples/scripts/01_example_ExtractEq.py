# %% Run the function to extract the data from local files
import os
from obspy import read
from os.path import exists
from os import mkdir
from rfsed.ExtractEq import ExtractEq
staname='NE301'
stalat=6.183432
stalon=53.48334
current_dir = os.path.dirname(os.path.realpath(__file__))
target_dir2 = os.path.sep.join(current_dir.split(os.path.sep)[:-2])
datapath= target_dir2 + '/data/Raw_Data_DeepNL_Groningen/' 
target_dir1=os.path.sep.join(current_dir.split(os.path.sep)[:-1])
savepath = target_dir1 + '/earthquake_data/Extracted_Eq_Data/'
if not os.path.exists(savepath):
    os.mkdir(savepath)
filename = savepath + 'eq_data' + '.h5'
catalog = target_dir1 + "/catalog/eq_events.quakeml"
Request_window=[-50, 150] #Time relative to first P arrival
ExtractData= ExtractEq(datapath, filename, catalog, stalat, stalon, Request_window)
#%% Check the extracted data
ExtractedEQ_Data=read(filename )
print(len(ExtractedEQ_Data), ' traces in the extracted data')
print(ExtractedEQ_Data)
# %%
