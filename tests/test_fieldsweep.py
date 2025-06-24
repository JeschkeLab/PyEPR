
from pyepr.fieldsweep_analysis import *
from pyepr.dataset import *
import pytest

def test_create_Nmodel():
    model = create_Nmodel(34*1e3)

    assert model is not None
    assert hasattr(model, "Boffset")
    assert hasattr(model, "gy")
    assert hasattr(model, "gz")
    assert hasattr(model, "GB")
    assert hasattr(model, "az")
    assert hasattr(model, "axy")
    assert hasattr(model, "scale")

def test_create_Nmodel_with_params():
    model = create_Nmodel(34*1e3)
    B_axis = np.linspace(12100, 12300, 10)
    params = {"B":B_axis/1e1,"az":3.66,"axy":0.488,"gy":2.01,"gz":2.0,"GB":0.45,"scale":1, "Boffset":0.7}
    V = model(**params)
    ref_values = np.array([1.68557056e-11, 4.50164637e-01, 8.92010015e-01, 7.09120791e-01,
       4.82187836e-01, 3.09046641e-01, 1.07465055e-01, 1.33308381e-02,
       3.41919276e-13, 2.24934385e-34])
    assert np.allclose(V,ref_values) # refvalues from 15th Nov 2023
