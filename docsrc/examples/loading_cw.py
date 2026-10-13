"""
Loading a CW EPR file
===============================================

PyEPR can also be used to load CW EPR data files, for plotting and analysis in python.

At the moment, PyEPR is focused on providing analysis for pulse EPR data, and the CW EPR functionality is limited to loading and plotting data.

Supported file formats include:
- Bruker .DTA/.DSC files
- Magnetech (binary) .dat files
- PyEPR .h5 files (saved using the ``dataset.epr.save()`` method)

"""

import pyepr as epr
import matplotlib.pyplot as plt

