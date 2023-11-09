from rfsed.rf.rfstream import read_rf, rfstats
from pkg_resources import resource_filename
import numpy as np
import os


def rfMoho_example():
    """
    Return low-frequency receiver functions in a stream by read_rf().
    """
    fname = resource_filename('rfsed', 'data/rfstreams_Moho/rfstreams.h5')
    stream = read_rf(fname)  
    return stream

def rfSed_example():
    """
    Return high-frequency receiver functions in a stream read_rf().
    """
    fname = resource_filename('rfsed', 'data/rfstreams_Sed/rfstreams.h5')
    stream = read_rf(fname)  
    return stream

def raw_data_example():
    """
    Return directory to raw seismic record data.
    """
    data_dir = resource_filename('rfsed', 'data/Raw_Data_DeepNL_Groningen/')

    return data_dir

def catalog_example():
    """
    Return directory to earthquake catalog.
    """
    catalog_dir = resource_filename('rfsed', 'data/Earthquake_Catalog/eq_events.quakeml')
    return catalog_dir

def save_Eq_data():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = resource_filename('rfsed', 'examples/earthquake_data/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_catalog():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = resource_filename('rfsed', 'examples/catalog/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_plot():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = resource_filename('rfsed', 'examples/plots/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_tests():
    """
    Return directory to save test plots.
    """
    savepath = resource_filename('rfsed', 'tests/plot/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath