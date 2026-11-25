from pyepr.pulses import RectPulse, GaussianPulse, Detection
import pytest 
import numpy as np
import pyepr as epr


def test_gaussian_pulse():
    pulse = GaussianPulse(tp=100, freq=0, scale=0.5)

    pulse.func()