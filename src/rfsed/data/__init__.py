# Copyright (c) 2023 Stephen Akinremi MIT license
"""
rfsed Documentation
================

rfsed is a Python framework for analysizing receiver functions.
It is developed specifically to implement different techniques of analysing receiver functions
from stations overlying sedimentary layers. Receiver functions techniques implemented in rfsed are:
1. H-k stacking of Zhu and Kanamori (2000)
2. Sequential H-k stacking of Yeck et al., (2013)
3. Resonance filtering and modified H-k stacking of Yu et al., (2015)
4. H-k stacking and Waveform Fitting two-step method of Akinremi et al., (2024)
6. Analysis of the synthetic reciever functions with the above-mentioned methods.

Beside these receiver function methods, rfsed has the following features:
1. Creating publication quality figures for the results of the analysis.
2. Extracting earthquake data from local seismic record files.
3. Multiprocessing options for waveform fitting and extracting earthquake data from local seismic record files.


Receiver functions streams are handled by the 'RFStream' class of the rf open software, and it inherits a lot of useful methods from the Obspy class 'Stream'.
It is supported via the obspyh5 package. For more information on class "RFStream", see documentation on rf (https://rf.readthedocs.io/en/latest/). 
The rf framework is included in the rfsed package for completness under The MIT License (MIT) Copyright (c) 2013-2019 Tom Eulenfeld.
In the module to extract earthquake data from local seismic records, read and write support for necessary metadata is provided 
for SAC, SeismicHanlder and HDF5 formats based on ObsPy.

The forward receiver function calculation is done using SynR, which is modified after the SynSeis module of Mijian Xu's seispy project (under the GNU GENERAL PUBLIC LICENSE).

Method
------
The receiver function technique is a well-established technique to image velocity contrasts in the subsurface by using isolated P–to-S wave 
conversions (P-receiver functions) and the reverberations generated at such discontinuities. The velocity contrast at the crust-mantle 
boundary (Moho) represents a major discontinuity that receiver functions are highly sensitive to, making them useful to constrain average 
crustal thickness and Vp/Vs. Several methods such as H-k stacking, has been developed to investigate the crustal thickness and Vp/Vs ratio 
using receiver functions. However, the presence of near-surface low-velocity sedimentary layer can obscure Moho signals due to the additional 
P-to-S wave conversions and  associated reverberated phases created at the sediment-basement discontinuity. These additional intra-crustal phases 
can have large amplitudes and similar arrival times to Moho signals, which makes it difficult to retrieve Moho information using standard RF techniques.
Several other methods have been developed to deal with the complex RF in sedimentary basins and recover crustal strucures (crustal thickness and Vp/Vs).
These methods include: Iterative trial and error fitting of observed receiver functions with synthetic ones,  wavefield continuation and decomposition technique, 
application of band-pass or resonance filters to remove the sediment signals in the receiver functions, sequential H-k stacking, and combination of H-k stacking 
and waveform fitting. 


Installation
------------

Dependencies of rfsed are
    * ObsPy_ and some of its dependencies,
    * cartopy, geographiclib, shapely,
    * mtspec_ for multitaper deconvolution,
    * toeplitz_ for faster time domain deconvolution (optional),
    * obspyh5_ for hdf5 file support,
    * tqdm for progress bar in batch processing (optional).

The easiest way to install rfsed is via `pip`_::
    
        pip install rfsed

Here are some instructions to install rfsed into a fresh conda environment::

    conda create -n rfsed python==3.9 
    conda activate rfsed
    git clone https://github.com/akinremisa/rfsed.git
    cd rfsed #Change directory to the same directory that this repo is in
    pip install -e .

Alternatively, the environment can be created and the package installed using the environment.yml file in the root directory of this repo::

1. Clone the rffw repository
git clone https://github.com/akinremisa/rfsed.git
2. Change directory to the same directory that this repo is cloned (i.e., same directory as environment.yml)
cd rfsed
3. Create the evniroment and install the rfsed package using the environment.yml file
conda env create -f environment.yml
4. Activate the environement
conda activate rfsed


Using the Python module/Tutorial
--------------------------------
Examples of how to use the methods in this Python module can be found in the `examples` directory.

Miscellaneous
-------------

Please feel free to request features, report bugs or contribute code on
GitHub_. The code is continuously tested by Github Actions.

Citation
--------

If you found this package useful, please consider citing it.
.....................................


.. _ObsPy: http://www.obspy.org/
.. _pip: http://www.pip-installer.org/
.. _obspyh5: https://github.com/trichter/obspyh5/
.. _toeplitz: https://github.com/trichter/toeplitz/
.. _mtspec: https://github.com/krischer/mtspec


.. _GitHub: https://github.com/akinremisa/rfsed/
"""
