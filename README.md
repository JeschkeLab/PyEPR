# PyEPR

![PyPI - Version](https://img.shields.io/pypi/v/pyEPR-ESR)
[![DOI](https://zenodo.org/badge/888368760.svg)](https://doi.org/10.5281/zenodo.17107010)



## About
PyEPR is a Python based software package for designing, running and processing EPR experiments on a wide variety of hardware. It is built to be modular and extensible, allowing easy integration of new hardware and experiment types. PyEPR is developed by the Jeschke Lab @ ETH. To extend the functionality of PyEPR we are actively looking for contributors.

PyEPR currently supports BRUKER ElexSys spectrometers running Xepr 2.9 and Andrin Doll style homebuilt spectrometers running Matlab based control software. Extending support for new hardware is straightforward and we are happy to assist anyone interested in doing so.

PyEPR was developed to provide the hardware abstraction layer for [autoDEER](https://jeschkelab.github.io/autoDEER/), a software package for running automated DEER experiments. However, PyEPR can also be used as a general purpose EPR data processing and simulation library.

## Features
- Fully python based, open-source and free to use
- Intuitive object-oriented pulse sequencer
- Pre-defined common EPR experiments (CW, Hahn Echo, Inversion Recovery, Carr-Purcell, DEER, etc.)
- Easy to define custom experiments
- Pre-defined common pulse shapes (rectangular, Gaussian, sech/tanh, etc.)
- Easy to define custom pulse shapes
- Hardware abstraction layer for interfacing with different spectrometers
- BRUKER PulseSpel compiler from PyEPR sequences

## Installation
It is recommended to install PyEPR in a virtual environment, from source.

## Dependencies
- numpy
- scipy
- matplotlib
- [DeerLab](https://jeschkelab.github.io/DeerLab/)
- pyyaml
- xarray
- h5netcdf
- toml
- numba

## Requirements
PyEPR is generally compatible with Windows, Mac and Linux and requires Python 3.11 or 3.12. 
Support for new Python versions will only be added once they become stable, and have been fully tested with PyEPR. This can take some time, so please be patient.
The specific hardware implementation may add additional limitations for example, when using XeprAPI a modern Linux OS is required.

## Citing PyEPR
At the moment, no paper is associated with PyEPR. Once published, it will be linked here. It is kindly requested that in the meantime the appropriate DOI number is cited. The associated Zotero entry is available [here](missing).

Citing academic software is important as it helps to ensure the long-term sustainability of the software, and allows the developers to track the impact of their work and secure future funding. It also helps to provide credit to the developers for their hard work.

## Contributing
Contributions to PyEPR are welcome! If you have discovered an issue or have a feature request, please open an issue on the GitHub repository. If you would like to contribute code, please fork the repository and submit a pull request. If you have any questions or need help, please open an issue or contact the authors.

## License
PyEPR is licensed under the GNU GPLv3 public license, and is released without
warranty or liability. Commercial use is allowed, however it is advised to contact the authors for support.

Copyright © 2021-2025: Hugo Karas, Stefan Stoll and Gunnar Jeschke
