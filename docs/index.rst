
.. figure:: ./logo/rfsed_logo.png
   :align: center

Documentation
=============

rfsed is developed specifically to implement different techniques of **analysing receiver functions
from stations overlying sedimentary layers.** The package is **adaptable, efficient, and easy-to-use** for different analysis receiver functions
obtained from stations overlying sedimentary layers.

Receiver functions techniques implemented in rfsed are:
+ H-k stacking of Zhu and Kanamori (2000)
+ Sequential H-k stacking of Yeck et al., (2013)
+ Resonance filtering and modified H-k stacking of Yu et al., (2015)
+ H-k stacking and Waveform Fitting two-step method of Akinremi et al., (2024)
+ Analysis of the synthetic reciever functions with the above-mentioned methods.

Beside these receiver function methods, rfsed has the following features:
+ Creating publication quality figures for the results of the analysis.
+ Extracting earthquake data from local seismic record files.
+ Multiprocessing options for waveform fitting and extracting earthquake data from local seismic record files.


Receiver functions streams in the rfsed are handled by the 'RFStream' class of the [rf](https://github.com/trichter/rf) open software, and it inherits a lot of useful methods from the Obspy class 'Stream'. It is supported via the obspyh5 package. For more information on class "RFStream", see documentation on [rf](https://rf.readthedocs.io/en/latest/). 
In the rfsed modules to extract earthquake data from local seismic records, read and write support for necessary metadata is provided for SAC, SeismicHanlder and HDF5 formats based on [ObsPy](https://github.com/obspy/obspy).

This package is fully in ``python`` and common processing workflows are described in the ``Jupyter`` notebooks accompanying this package.

.. toctree::
   :maxdepth: 1
   :caption: Quick Links

   links

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   api

.. toctree::
   :maxdepth: 1
   :caption: Module

   elast
   utils
   wiggles

.. toctree::
   :maxdepth: 1
   :caption: Jupyter Notebooks

   Example 1: A simple example for OBS <https://nbviewer.jupyter.org/github/paudetseis/Telewavesim/blob/master/telewavesim/examples/Notebooks/sim_obs_Audet2016.ipynb>
   Example 2: Lower crustal anisotropy <https://nbviewer.jupyter.org/github/paudetseis/Telewavesim/blob/master/telewavesim/examples/Notebooks/sim_Prfs_Porter2011.ipynb>
   Example 3: Simulating SKS splitting <https://nbviewer.jupyter.org/github/paudetseis/Telewavesim/blob/master/telewavesim/examples/Notebooks/sim_SKS.ipynb>

