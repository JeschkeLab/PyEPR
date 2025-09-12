
Welcome to PyEPR!
====================================

.. caution:: 
    This documentation is a work in progress and not yet complete. 
    Please be patient and if anything is unclear and please feel free to contact the developers.

PyEPR is a Python package for designing for building and running automated Electron Paramagnetic Resonance (EPR) sequences.
The package has been designed to support a wide range of EPR spectrometers, both commercial and home-built.


PyEPR's Key Features
-----------------------

- Fully python based, open-source and free to use
- Intuitive object-oriented pulse sequencer
- Pre-defined common EPR experiments (CW, Hahn Echo, Inversion Recovery, Carr-Purcell, DEER, etc.)
- Easy to define custom experiments
- Pre-defined common pulse shapes (rectangular, Gaussian, sech/tanh, etc.)
- Easy to define custom pulse shapes
- Hardware abstraction layer for interfacing with different spectrometers
- BRUKER PulseSpel compiler from PyEPR sequences

.. warning:: 
    PyEPR is an actively developed software package, that is still very much a work in process. Please consider this to be a beta release.
      
.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: User Guide

    ./install.rst
    ./tutorial.rst
    ./API_docs.rst

.. toctree::
    :hidden:
    :caption: About

    ./releasenotes.rst
    ./contributing.rst
    Github <https://github.com/JeschkeLab/PyEPR>
    autoDEER <https://github.com/JeschkeLab/autoDEER>


