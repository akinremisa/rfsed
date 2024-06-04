---
title: 'rfsed: Receiver function analysis and dealing with sediment effects'
tags:
  - Python
  - geophysics
  - seismology
  - receiver functions
  - sediment effects
  - Moho discontinuity
authors:
  - name: Stephen Akinremi
    orcid: 0000-0001-9162-7156
    affiliation: 1
  - name: Islam Fadel
    orcid: 0000-0002-0091-8175
    affiliation: 1
  - name: Mark van der Meijde
    orcid: 0000-0002-8762-585X
    affiliation: 1
affiliations:
 - name: Faculty of Geo-Information Science and Earth Observation (ITC), University of Twente, Enschede, The Netherlands
   index: 1
date: 2024
bibliography: paper.bib

---

# Introduction

The receiver function technique is a well-established technique to image velocity contrasts in the subsurface (such as the crust-mantle boundary - Moho) using isolated P–to-S wave conversions (P-receiver functions) and the reverberations generated at such discontinuities. Several methods such as H-k stacking [@Zhu2000], have been developed to investigate the average crustal thickness and Vp/Vs ratio using receiver functions. However, the presence of a near-surface low-velocity sedimentary layer can obscure Moho phases due to the additional P-to-S wave conversions and associated reverberated phases created at the sediment-basement discontinuity. These additional intra-crustal phases can have large amplitudes and similar arrival times as the Moho phases, which makes it difficult to retrieve Moho information using standard receiver function techniques. 

# Statement of need
``rfsed`` is a Python software for receiver function analysis that includes methods for dealing with sediment effects. ``rfsed`` presents a new approach for retrieving reliable crustal thickness and Vp/Vs from stations overlying sedimentary layer. The techniques derives sediment thickness and Vp/Vs using H-K stacking of the high-frequency receiver function, followed by a waveform fitting approach to retrieve the average crustal thickness and Vp/Vs [@Akinremi2024]. Moreover, ``rfsed`` contains implementations of the most common receiver function approaches for dealing with the sediment effect with possible synthetic testing capabilities. ``rfsed`` has already been used in [@Akinremi2024] for multimethods imaging of the crustal structure using receiver functions obtained from stations overlying sedimentary layer.

# Key Functionality
``rfsed`` contains modules to carry out H-k stacking [@Zhu2000], sequential H-k stacking [@Yeck2013], resonance filtering and modified H-k stacking [@Yu2015], and waveform fitting [@Akinremi2024] with possible synthetic waveform generation for 1D earth models to test the different methods. It comes with tools to create high-quality figures, which include result plots for H-k stacking, sequential H-k stacking (e.g., Figure 1), resonance filtering, and waveform fitting methods (e.g., Figure 2). Besides these methods, ``rfsed`` has modules for extracting earthquake waveforms from local seismic record files. There are multiprocessing options for waveform fitting and extracting earthquake data from local seismic record files for higher efficiency.

![Example of a sequential H-k stacking plot for receiver functions obtained from station ROLD (Network: NL) (a) sediment layer (b) Moho layer, generated using ``rfsed``](paper_figures/Figure1.png)

``rfsed`` is adaptable, efficient, and easy-to-use by both researchers and students. Receiver function streams in ``rfsed`` are handled by the 'RFStream' class of ``rf`` [@Eulenfeld2020]. ``rfsed`` can be installed from [PyPI](https://pypi.org/project/rfsed/). Online documentation and tutorials are available on the [project site](https://akinremisa.github.io/rfsed/).

![Example of a waveform-fitting result plot generated using ``rfsed``](paper_figures/Figure2.png)

# Availability

The software is distributed under a BSD License and is available from [rfsed](https://github.com/akinremisa/rfsed).


# Acknowledgements

``rfsed`` was initiated as part of the DeepNL-DICTUM project supported by the Nederlandse Organisatie voor Wetenschappelijk Onderzoek (NWO) DeepNL program grant number DEEP.NL.2020.010. The data used in creating Figure 1 is from the Netherlands Seismic and Acoustic Network operated by Royal Netherlands Meteorological Institute (KNMI), with network code NL. 



# References
