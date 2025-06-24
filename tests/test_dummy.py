from pyepr.hardware.dummy import *
import pytest

def test_create_dummyInterface():
    config = {}
    config['Spectrometer'] = {'Dummy': {'speedup': 100, 'SNR': 1, 'ESEEM_depth': 1, 'noise_level': 0.001, 'Sample': {'name': 'test', 'conc': 1}}, 'Bridge': {}}
    config['Resonators'] = {'test': {'Center Freq': 34.0, 'Q': 100, 'nu1': 75}}
    dummy = dummyInterface(config)
    
    assert dummy is not None

    pytest.approx(dummy.mode(34.0),75)