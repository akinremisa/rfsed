.. figure:: ./logo/rfsed_logo_horizontal.png
   :align: center

Documentation
================

**rfsed** is developed specifically to implement different techniques of 
**analysing receiver functions from stations overlying sedimentary layers.** 
The software is **adaptable, efficient, and easy-to-use** for different 
analysis of receiver functions obtained from stations overlying sedimentary layers.


Receiver functions techniques implemented in **rfsed** are:

* H-k stacking (one layer) of Zhu and Kanamori (2000)

* Sequential H-k stacking (two layers) of Yeck et al., (2013)

* Resonance filtering and modified H-k stacking of Yu et al., (2015)

* H-k stacking and Waveform Fitting (two-step method) of Akinremi et al., (2024)

* Analysis of the synthetic reciever functions with the above-mentioned methods.


Beside these receiver function methods, rfsed has the following features:

* Extracting earthquake data from local seismic record files.

* Multiprocessing options for waveform fitting and extracting earthquake data from local seismic record files.

* Creating publication quality figures for the results of the analysis.


Receiver functions streams in the rfsed are handled by the `RFStream` class of the 
`rf <https://github.com/trichter/rf>`_ open software, and it inherits a lot of useful 
methods from the Obspy class `Stream`. It is supported via the obspyh5 package. 
For more information on class `RFStream`, see documentation on 
`rf <https://github.com/trichter/rf>`_. In the rfsed modules to extract earthquake 
data from local seismic records, read and write support for necessary metadata is 
provided for SAC, SeismicHanlder and HDF5 formats based on 
`ObsPy <https://github.com/obspy/obspy>`_.

This package is fully in ``python`` and common processing workflows are described in 
the ``Jupyter`` notebooks accompanying this package.


.. toctree::
   :maxdepth: 1
   :caption: Quick Links

   links

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   gettingstarted
   
.. toctree::
   :maxdepth: 1
   :caption: Module

   hkZhu
   hkSeqYeck
   ReverbFilter
   hkYu
   synrf
   WaveformFitting
   WaveformFittingMultiproc
   ExtractEq
   ExtractEqMultiproc
   

   .. toctree::
   :maxdepth: 1
   :caption: Jupyter Notebooks


   `Example 1: Downloading Earthquake Waveform Data from IRIS <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example1_Download_Waveform_Data_from_IRIS.ipynb>`_
   `Example 2: Calculation of Receiver Functions <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example2_Calculate_Receiver_Functions.ipynb>`_
   `Example 3: Download Catalog of earthquake events <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example3_example_Get_Catalog.ipynb>`_
   `Example 4: Example of Extracting Earthquake Waveform Data from seismic data from local drive <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example4_Extract_Earthquake_Waveform_from_Local_Data.ipynb>`_
   `Example 5: Example of Extracting Earthquake Waveform Data from seismic data from local drive - Using Multiple Processors <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example5_Extract_Earthquake_Waveform%28Multiprocessing%29_from_Local_Data.ipynb>`_
   `Example 6: H-K Stacking of Receiver Functions <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example6_One_layer_H-K_Stacking_hkZhu.ipynb>`_
   `Example 7: Sequential H-K Stacking of Receiver Functions <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example7_Sequential_H-K_Stacking_hkYeck_Sequential.ipynb>`_
   `Example 8: Resonance Filtering and Modified H-K Stacking of Receiver Functions <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example8_Resonance_Filtering_and_Modified_H-K_Stacking_hkYu.ipynb>`_
   `Example 9: Calculating Synthetic Reciever Function <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example9_Calculating_Synthetic_Receiver_Functions.ipynb>`_
   `Example 10: Calculation and Analysis of Synthetic Receiver Functions <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example10_Synthethic_Receiver_Function_Analysis.ipynb>`_
   `Example 11: Waveform Fitting of Receiver Functions <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example11_Waveform_Fitting_of_Receiver_Functions.ipynb>`_
   `Example 12: Waveform Fitting of Receiver Functions - Using Multiple Processors <https://nbviewer.org/github/akinremisa/rfsed/blob/main/src/rfsed/examples/notebooks/Example12_Waveform_Fitting%28Multiprocessing%29_of_Receiver_Functions.ipynb.ipynb>`_
 