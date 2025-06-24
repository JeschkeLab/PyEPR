import pytest
import numpy as np


from pyepr.resonator_profile_analysis import *

def test_resonator_profile_analysis__dummy_create():
    fc = 34.0
    Q = 100
    nu_max = 50

    freq_axis = np.linspace(fc-0.2,fc+0.2,20)

    analysis = ResonatorProfileAnalysis._dummy_create(freq_axis,nu_max,Q,fc)