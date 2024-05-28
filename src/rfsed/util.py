from rf.rfstream import read_rf, rfstats
from pkg_resources import resource_filename
import numpy as np
import os
import pkg_resources as _pkg_resources
from distutils import dir_util as _dir_util


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

def read_raw_waveform_data():
    """
    Return streams of raw waeform data downloaded from IRIS using read_rf().
    """
    fname = resource_filename('rfsed', 'examples/IRIS_Downloaded_Waveform/00_eq_data.h5')
    stream = read_rf(fname)  
    return stream

def save_Eq_data():
    """
    Return directory to save extracted earthquake data.
    """
    savepath = resource_filename('rfsed', 'examples/earthquake_data/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_IRIS_waveform():
    """
    Return directory to save downloaded earthquake waveform data.
    """
    savepath = resource_filename('rfsed', 'examples/IRIS_Downloaded_Waveform/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_calculated_RF():
    """
    Return directory to save calculated receiver function.
    """
    savepath = resource_filename('rfsed', 'examples/Calculated_Receiver_Functions/')
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
    savepath = resource_filename('rfsed', 'examples/test_output/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def install_examples(path="./rfsed_examples"):
    """
    Install the examples notebooks and data files for rfsed in the given location.

    WARNING: If the path exists, the files will be written into the path
    and will overwrite any existing files with which they collide.

    """
    subpath_example = os.path.join(path, "examples")

    examples_path = _pkg_resources.resource_filename(
        "rfsed", "examples")

    _dir_util.copy_tree(
        examples_path,
        subpath_example,
        preserve_mode=1,
        preserve_times=1,
        preserve_symlinks=1,
        update=0,
        verbose=1,
        dry_run=0)
    
    subpath_data = os.path.join(path, "data")
    data_path = _pkg_resources.resource_filename(
        "rfsed", "data")

    _dir_util.copy_tree(
        data_path,
        subpath_data,
        preserve_mode=1,
        preserve_times=1,
        preserve_symlinks=1,
        update=0,
        verbose=1,
        dry_run=0)
