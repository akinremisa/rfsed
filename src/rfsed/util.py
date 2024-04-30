from rf.rfstream import read_rf, rfstats
import numpy as np
import os
current_directory = os.path.dirname(os.path.abspath(__file__))
main_dir = os.path.abspath(os.path.join(current_directory, '..', '..'))

def rfMoho_example():
    """
    Return low-frequency receiver functions in a stream by read_rf().
    """
    fname = main_dir + '/data/rfstreams_Moho/rfstreams.h5'
    stream = read_rf(fname)  
    return stream

def rfSed_example():
    """
    Return high-frequency receiver functions in a stream read_rf().
    """
    fname = main_dir + '/data/rfstreams_Sed/rfstreams.h5'
    stream = read_rf(fname)  
    return stream

def raw_data_example():
    """
    Return directory to raw seismic record data.
    """
    data_dir = main_dir + '/data/Raw_Data_DeepNL_Groningen/'

    return data_dir

def catalog_example():
    """
    Return directory to earthquake catalog.
    """
    catalog_dir = main_dir + '/data/Earthquake_Catalog/eq_events.quakeml'
    return catalog_dir

def read_raw_waveform_data():
    """
    Return streams of raw waeform data downloaded from IRIS using read_rf().
    """
    fname = main_dir + '/examples/IRIS_Downloaded_Waveform/00_eq_data.h5'
    stream = read_rf(fname)  
    return stream

def save_Eq_data():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = main_dir +  '/examples/earthquake_data/'
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_IRIS_waveform():
    """
    Return directory to save downloaded earthquake waveform data.
    """
    savepath = main_dir + '/examples/IRIS_Downloaded_Waveform/'
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_calculated_RF():
    """
    Return directory to save calculated receiver function.
    """
    savepath = main_dir + '/examples/Calculated_Receiver_Functions/'
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_catalog():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = main_dir + '/examples/catalog/'
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_plot():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = main_dir + '/examples/plots/'
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_tests():
    """
    Return directory to save test plots.
    """
    savepath = main_dir + '/tests/plots/'
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath