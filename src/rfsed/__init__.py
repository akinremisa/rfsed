# Copyright (c) 2023 Stephen Akinremi, Islam Fadel, MIT license
"""
rfsed Documentation
================

rfsed is a Python framework for analysizing receiver functions.
It is developed specifically to implement different techniques of analysing receiver functions
from stations overlying sedimentary layers. Receiver functions techniques implemented in rfsed are:
1. H-k stacking of Zhu and Kanamori (2000)
2. Sequential H-k stacking of Yeck et al., (2013)
3. Resonance filtering and modified H-k stacking of Yu et al., (2015)
4. Waveform Fitting grid search
5. Calculation of synthetic receiver functions using hrftn96 codes of Herrmann (2013)
6. Analysis of the synthetic reciever functions with the above mentioned methods.

Beside these receiver function methods, rfsed has the following features:
1. Creating publication quality figures for the results of the analysis.
2. Extracting earthquake data from local seismic record files.
3. Parallel process for waveform fitting grid search and extracting earthquake data from local seismic record files.


Receiver functions streams are handled by the 'RFStream' class of the rf open software, and it inherits a lot of useful methods from the Obspy class 'Stream'.
It is supported via the obspyh5 package. For more information on class "RFStream", see documentation on rf (https://rf.readthedocs.io/en/latest/). 
The rf framework is included in the rfsed package for completness under The MIT License (MIT) Copyright (c) 2013-2019 Tom Eulenfeld.
In the module to extract earthquake data from local seismic records, read and write support for necessary metadata is provided 
for SAC, SeismicHanlder and HDF5 formats based on ObsPy.

Method
------

...................
...................

Installation
------------

Dependencies of rf are

    * ObsPy_ and some of its dependencies,
    * obspyh5
    * cartopy, geographiclib, shapely,
    * mtspec_ for multitaper deconvolution,
    * toeplitz_ for faster time domain deconvolution (optional),
    * obspyh5_ for hdf5 file support (optional),
    * tqdm for progress bar in batch processing (optional).




Here are some instructions to install rf into a fresh conda environment::

    conda create -n rfsed python==3.9 
    conda activate rfsed
    git clone https://github.com/akinremisa/rfsed.git
    cd rfsed #Change directory to the same directory that this repo is in
    pip install -e .
   



Miscellaneous
-------------

Please feel free to request features, report bugs or contribute code on
GitHub_. The code is continuously tested by Github Actions.

Citation
--------

If you found this package useful, please consider citing it.
.....................................

.. _this: https://doi.org/10.17169/refubium-14424
.. _ObsPy: http://www.obspy.org/
.. _pip: http://www.pip-installer.org/
.. _obspyh5: https://github.com/trichter/obspyh5/
.. _toeplitz: https://github.com/trichter/toeplitz/
.. _mtspec: https://github.com/krischer/mtspec


.. _GitHub: https://github.com/akinremisa/rfsed/
"""
