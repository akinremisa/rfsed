<img src="https://github.com/akinremisa/rfsed/blob/main/doc/logo/rfsed_logo_horizontal.png" alt="rfsed logo" width="600"/>

## Receiver functions analysis and dealing with sediment effects in receiver functions.
[![python version](https://img.shields.io/pypi/pyversions/rf.svg)](https://python.org)




## A python package for receiver functions analysis and dealing with sediment effects 

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
The rf framework is included in the rfsed package for completness under The MIT License (MIT) Copyright (c) 2013-2019 Tom Eulenfeld. In the rfsed modules to extract earthquake data from local seismic records, read and write support for necessary metadata is provided for SAC, SeismicHanlder and HDF5 formats based on [ObsPy](https://github.com/obspy/obspy).

The receiver function forward calculation in rfsed is done using SynR, which is modified after the [SynSeis module - seispy](https://github.com/xumi1993/seispy) project (under the GNU GENERAL PUBLIC LICENSE).




## Installation and testing of this package

### Installation from PyPi
The easiest way to install rfsed is via `pip`_::

```bash
pip install rfsed
```

### Installation from source code using environment.yml
To obtain the latest updates, you can install rfsed from the source code, available on GitHub.

```bash
# Clone the rfsed repository from GitHub
git clone https://github.com/akinremisa/rfsed.git

# Change directory to the same directory that this repo is in (i.e., same directory as setup.py and environment.yml)
cd cd rfsed  # That's the standard name the folder should have

# Create the conda environment, install dependencies, and install the rfsed package using environment.yml
conda env create -f environment.yml

# Activate the conda environment
conda activate rfsed

```
### Test the rfsed package
You can test the package using pytest by running this command in the directory that has the '/tests/' folder

```bash
pytest -p no:logging tests
```
Or run individual tests in the '/tests/ directory

## Getting started
Access rfsed's documentation [here](https://rfsed.github.io/rfsed/).

rfsed comes with tutorials that demonstrates all its methods. You can find those in the `examples/` directory.

## Reporting Bugs / Contact the developers
This version is an early release of rfsed. If you encounter any issues or unexpected behaviour, please [open an issue](https://github.com/akinremisa/rfsed/issues/new) on GitHub.

## Questions?
If you have any questions about the package, please use the [discussions feature](https://github.com/akinremisa/rfsed/discussions/new/choose)

## Contributing
All contributions are welcome ... e.g. report bugs, discuss or add new features.

## Citing rfsed
If you found this package useful, please consider citing it.

#### References
+ Yeck, W. L., Sheehan, A. F., & Schulte-Pelkum, V. (2013). Sequential h-k stacking to obtain accurate crustal thicknesses beneath sedimentary basins. Bulletin of the Seismological Society of America, 103 , 2142-2150. doi: 10.1785/0120120290
+ Yu, Y., Song, J., Liu, K. H., & Gao, S. S. (2015). Determining crustal structure beneath seismic stations overlying a low-velocity sedimentary layer using receiver functions. Journal of Geophysical Research: Solid Earth, 120 , 3208-3218. doi: 10.1002/2014JB011610
+  Zhu, L., & Kanamori, H. (2000). Moho depth variation in southern California from teleseismic receiver functions. Journal of Geophysical Research: Solid Earth, 105, 2969-2980. doi: 10.1029/1999jb900322

##### Related receiver function projects
+ [rf](https://github.com/trichter/rf) including calculation of receiver functions
+ [seispy](https://github.com/xumi1993/seispy) including hk-stacking
+ [RFPy](https://github.com/paudetseis/RfPy) including hk-stacking, harmonic decomposition
+ [BayHunter](https://github.com/jenndrei/BayHunter) inversion of receiver functions and surface wave dispersion
+ [telewavesim](https://github.com/paudetseis/Telewavesim) synthetics
+ [PyGLImER](https://github.com/PyGLImER/PyGLImER) including common conversion point imaging
