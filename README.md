<img src="docs/logo/rfsed_logo_horizontal.png" alt="rfsed logo" width="600"/>


## A software for receiver functions analysis and dealing with sediment effects 

[![python version](https://img.shields.io/pypi/pyversions/rf.svg)](https://python.org)


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


Receiver functions streams in the **rfsed** are handled by the 'RFStream' class of the 
[rf](https://github.com/trichter/rf) open software, and it inherits a lot of useful 
methods from the Obspy class 'Stream'. It is supported via the obspyh5 package. 
For more information on class "RFStream", see documentation on 
[rf](https://rf.readthedocs.io/en/latest/). In the **rfsed** modules to extract earthquake 
data from local seismic records, read and write support for necessary metadata is 
provided for SAC, SeismicHanlder and HDF5 formats based on 
[ObsPy](https://github.com/obspy/obspy).




## Installation and testing of this software

### Installation from PyPi
The easiest way to install **rfsed** is via `pip`::

```bash
pip install rfsed
```

### Installation development version from source code
To obtain the latest updates, you can install **rfsed** from the source code from 
available on GitHub.

Clone the **rfsed** repository from GitHub
```bash
git clone https://github.com/akinremisa/rfsed.git
```
Change directory to the same directory that this repo is in
```bash
cd rfsed 
``` 
Installing using pip
```bash
pip install .
```
### Test the rfsed software
You can test the software using pytest by running this command in the software 
directory and will look for all available tests in the current directory and 
subdirectories recursively

```bash
pytest
```
Or run individual tests in the /tests/ directory

## Getting started
Access **rfsed**'s documentation [here](https://akinremisa.github.io/rfsed/).

**rfsed** comes with tutorials that demonstrates all its methods. See the documentation for  more details.

## Reporting Bugs / Contact the developers
This version is an early release of **rfsed**. If you encounter any issues or unexpected 
behaviour, please [open an issue](https://github.com/akinremisa/rfsed/issues/new) on GitHub.

## Questions?
If you have any questions about the software, please use the 
[discussions feature](https://github.com/akinremisa/rfsed/discussions/new/choose)

## Contributing
All contributions are welcome ... e.g. report bugs, discuss or add new features.

## Citing rfsed
If you use **rfsed** in your work, please consider citing the following paper.
+ Akinremi S., van der Meijde, M., Thomas, C., Afonso, J. C., Ruigrok, E., & Fadel, I. (2024). 
Waveform fitting of receiver functions for enhanced retrieval of crustal structure in the 
presence of sediments. Journal of Geophysical Research: Solid Earth, 5(129). https://doi.org/10.1029/2023JB028393

#### References
+ Akinremi, S., van der Meijde, M., Thomas, C., Afonso J. C., Ruigrok E., Fadel, I. (2024). 
Waveform fitting of receiver functions for enhanced retrieval of crustal structure in the 
presence of sediments. Journal of Geophysical Research: Solid Earth, 5(129). https://doi.org/10.1029/2023JB028393
+ Tom Eulenfeld T., (2020). rf: Receiver function calculation in seismology. Journal of Open Source Software, 5(48), 
1808, https://doi.org/10.21105/joss.01808
+ Yeck, W. L., Sheehan, A. F., & Schulte-Pelkum, V. (2013). Sequential h-k stacking to obtain accurate 
crustal thicknesses beneath sedimentary basins. Bulletin of the Seismological Society of America, 103, 
2142-2150. https://doi.org/10.1785/0120120290
+ Yu, Y., Song, J., Liu, K. H., & Gao, S. S. (2015). Determining crustal structure beneath seismic 
stations overlying a low-velocity sedimentary layer using receiver functions. Journal of Geophysical 
Research: Solid Earth, 120 , 3208-3218. https://doi.org/10.1002/2014JB011610
+ Zhu, L., & Kanamori, H. (2000). Moho depth variation in southern California from teleseismic 
receiver functions. Journal of Geophysical Research: Solid Earth, 105, 2969-2980. https://doi.org/10.1029/1999jb900322

##### Related receiver function projects
+ [rf](https://github.com/trichter/rf) including calculation of receiver functions
+ [seispy](https://github.com/xumi1993/seispy) including hk-stacking
+ [RFPy](https://github.com/paudetseis/RfPy) including hk-stacking, harmonic decomposition
+ [BayHunter](https://github.com/jenndrei/BayHunter) inversion of receiver functions and surface wave dispersion
+ [telewavesim](https://github.com/paudetseis/Telewavesim) teleseismic body wave modeling through stacks of anisotropic layers
+ [PyGLImER](https://github.com/PyGLImER/PyGLImER) including common conversion point imaging
