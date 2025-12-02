Installation
------------

PyEPR can either be installed via pip or from source. 

It is recommended that virtual environments are used to install PyEPR, one recommended way to do this is with poetry.

To install PyEPR via pip, simply run the following command:

.. code-block:: bash

    pip install pyepr-esr

or with poetry:

.. code-block:: bash

    poetry add pyepr-esr

To install PyEPR from source, clone the repository and run the following command:

.. code-block:: bash
    
    git clone https://github.com/JeschkeLab/PyEPR

    pip install .

Requirements
++++++++++++

PyEPR requires:
    - Python >= 3.11 < 3.13
    - Numpy 
    - Scipy
    - Matplotlib
    - pyyaml
    - xarray
    - h5netcdf
    - toml
    - deerlab (https://github.com/JeschkeLab/DeerLab)
    - numba
    - psutil
