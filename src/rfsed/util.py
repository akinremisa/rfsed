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

from rf.rfstream import read_rf
from pkg_resources import resource_filename
import os
import pkg_resources as _pkg_resources
from distutils import dir_util as _dir_util

def rfMoho_example():
    """
    Return low-frequency receiver functions in a stream by read_rf().

    Example
    -------

    >>> # Initialize the rfMoho_example module:
    >>> from rfsed.util import rfMoho_example
    >>> rfMoho=rfMoho_example()
    """

    fname = resource_filename('rfsed', 'data/rfstreams_Moho/rfstreams.h5')
    stream = read_rf(fname)  
    return stream

def rfSed_example():
    """
    Return high-frequency receiver functions in a stream read_rf().
    
    Example
    -------

    >>> # Initialize the rfSed_example module:
    >>> from rfsed.util import rfSed_example
    >>> rfSed=rfSed_example()
    """

    fname = resource_filename('rfsed', 'data/rfstreams_Sed/rfstreams.h5')
    stream = read_rf(fname)  
    return stream

def raw_data_example():
    """
    Return directory to raw seismic record data.

    Example
    -------

    >>> # Initialize the raw_data_example module:
    >>> from rfsed.util import raw_data_example
    >>> data_dir=raw_data_example()
    """

    data_dir = resource_filename('rfsed', 'data/Raw_Data/')
  

    return data_dir

def catalog_example():
    """
    Return directory to earthquake catalog.

    Example
    -------

    >>> # Initialize the catalog_example module:
    >>> from rfsed.util import catalog_example
    >>> catalog_dir=catalog_example()
    """

    catalog_dir = resource_filename('rfsed', 'data/Earthquake_Catalog/eq_events.quakeml')
    return catalog_dir

def read_raw_waveform_data():
    """
    Return streams of raw waeform data downloaded from IRIS using read_rf().

    Example
    -------

    >>> # Initialize the read_raw_waveform_data module:
    >>> from rfsed.util import read_raw_waveform_data
    >>> rfstream=read_raw_waveform_data()
    """

    fname = resource_filename('rfsed', 'examples/IRIS_Downloaded_Waveform/00_eq_data.h5')
    stream = read_rf(fname)  
    return stream

def save_Eq_data():
    """
    Return directory to save extracted earthquake data.

    Example
    -------

    >>> # Initialize the save_Eq_data module:
    >>> from rfsed.util import save_Eq_data
    >>> savepath=save_Eq_data()
    """

    savepath = resource_filename('rfsed', 'examples/earthquake_data/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_IRIS_waveform():
    """
    Return directory to save downloaded earthquake waveform data.

    Example
    -------

    >>> # Initialize the save_IRIS_waveform module:
    >>> from rfsed.util import save_IRIS_waveform
    >>> savepath=save_IRIS_waveform()
    """

    savepath = resource_filename('rfsed', 'examples/IRIS_Downloaded_Waveform/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_calculated_RF():
    """
    Return directory to save calculated receiver function.

    Example
    -------

    >>> # Initialize the save_calculated_RF module:
    >>> from rfsed.util import save_calculated_RF
    >>> savepath=save_calculated_RF()
    """

    savepath = resource_filename('rfsed', 'examples/Calculated_Receiver_Functions/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_catalog():
    """
    Return directory to save extracted earthquake data.

    Example
    -------

    >>> # Initialize the save_catalog module:
    >>> from rfsed.util import save_catalog
    >>> savepath=save_catalog()
    """

    savepath = resource_filename('rfsed', 'examples/catalog/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_plot():
    """
    Return directory to save extracted earthquake data.

    Example
    -------

    >>> # Initialize the save_plot module:
    >>> from rfsed.util import save_plot
    >>> savepath=save_plot()
    """

    savepath = resource_filename('rfsed', 'examples/plots/')
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    return savepath

def save_tests():
    """
    Return directory to save test plots.

    Example
    -------

    >>> # Initialize the save_tests module:
    >>> from rfsed.util import save_tests
    >>> savepath=save_tests()
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

    :param path: The path to install the examples.
    :type path: str

    Example
    -------

    >>> # Initialize the install_examples module:
    >>> from rfsed.util import install_examples
    >>> install_examples(path='rfsed_examples')
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
