"""
Xepr API Interface
===============================================

Connecting to a Bruker Xepr (or Xenon) based spectrometer needs to be done using the XeprAPI.
Using this API, PyEPR can control the spectrometer and upload pulse sequences to it for execution.
It also has a lower-level interface which allows for more advanced control of the spectrometer, 
which can be useful for CW as well as pulse EPR experiments.

This code may work with Xepr > 2.7, but is only recommened for Xepr 2.9 and above.
"""

# %%
# Import packages including the Xepr interface
# --------------------------------------------

import pyepr as epr
# from pyepr.hardware.Bruker_AWG import BrukerAWG

