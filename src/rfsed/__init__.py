# Copyright (c) 2023 Stephen Akinremi BSD 3-Clause License
"""
rfsed Documentation
================
**rfsed** is developed specifically to implement different techniques of 
**analysing receiver functions from stations overlying sedimentary layers.** 
The software is **adaptable, efficient, and easy-to-use** for different 
analysis of receiver functions obtained from stations overlying sedimentary layers.

Receiver functions techniques implemented in **rfsed** are:
+ H-k stacking (one layer) of Zhu and Kanamori (2000)
+ Sequential H-k stacking (two layers) of Yeck et al., (2013)
+ Resonance filtering and modified H-k stacking of Yu et al., (2015)
+ H-k stacking and Waveform Fitting (two-step method) of Akinremi et al., (2024)
+ Analysis of the synthetic reciever functions with the above-mentioned methods.


Beside these receiver function methods, **rfsed** has the following features:
+ Extracting earthquake data from local seismic record files.
+ Multiprocessing options for waveform fitting and extracting earthquake 
data from local seismic record files.
+ Creating publication quality figures for the results of the analysis.


Receiver functions streams in the **rfsed** are handled by the `RFStream` class of the 
[rf](https://github.com/trichter/rf) open software, and it inherits a lot of useful 
methods from the Obspy class `Stream`. It is supported via the obspyh5 package. 
For more information on class `RFStream`, see documentation on 
[rf](https://rf.readthedocs.io/en/latest/). In the **rfsed** modules to extract earthquake 
data from local seismic records, read and write support for necessary metadata is 
provided for SAC, SeismicHanlder and HDF5 formats based on 
[ObsPy](https://github.com/obspy/obspy).
"""

from .__version__ import __version__