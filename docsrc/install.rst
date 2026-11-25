Installation
------------

PyEPR is published on PyPI as ``pyepr-esr`` and supports **Python 3.11 or newer**.
It can be installed from PyPI or from source.

PyEPR is a standard Python package, so any Python package manager can install it.
This page gives instructions for `uv <https://docs.astral.sh/uv/>`_, a fast
package and project manager that we use for development, and for plain ``pip``.
Other tools, such as Poetry, Pixi or conda/mamba (installing PyEPR with ``pip``
inside a conda environment), work too. Use whichever you are most comfortable with.

.. important::

    We recommend installing PyEPR in a **local virtual environment** with a
    **recent Python version** (3.11 or newer). Do not install it into the
    system Python that ships with your operating system: it is often outdated,
    and installing packages into it can conflict with system tools.

    A virtual environment keeps PyEPR and its dependencies separate from other
    projects. uv makes this easy because it can download a recent Python
    version for you (for example ``uv venv --python 3.12``), so you do not need
    to install Python yourself.

Installing uv (optional)
++++++++++++++++++++++++

You only need this step if you want to use uv. If you do not already have it,
install it with:

.. tab-set::
    :sync-group: os

    .. tab-item:: macOS / Linux
        :sync: unix

        .. code-block:: bash

            curl -LsSf https://astral.sh/uv/install.sh | sh

    .. tab-item:: Windows
        :sync: windows

        .. code-block:: powershell

            powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

For other install methods (Homebrew, pipx, WinGet, ...), see the
`uv installation guide <https://docs.astral.sh/uv/getting-started/installation/>`_.

Installing from PyPI
++++++++++++++++++++

.. tab-set::
    :sync-group: tool

    .. tab-item:: uv (project)
        :sync: uv

        To add PyEPR as a dependency of your own uv project (this creates the
        virtual environment automatically):

        .. code-block:: bash

            uv init my-epr-project   # skip if you already have a project
            cd my-epr-project
            uv add pyepr-esr

        Run scripts inside the project environment with ``uv run``:

        .. code-block:: bash

            uv run python my_script.py

    .. tab-item:: uv (pip interface)
        :sync: uv-pip

        To install into a standalone virtual environment:

        .. code-block:: bash

            uv venv --python 3.12
            source .venv/bin/activate      # Windows: .venv\Scripts\activate
            uv pip install pyepr-esr

    .. tab-item:: pip
        :sync: pip

        .. code-block:: bash

            python -m venv .venv
            source .venv/bin/activate      # Windows: .venv\Scripts\activate
            python -m pip install pyepr-esr

Installing from source
++++++++++++++++++++++

Install from source to get the latest development version, or if you want to
contribute to PyEPR. Clone the repository first:

.. code-block:: bash

    git clone https://github.com/JeschkeLab/PyEPR.git
    cd PyEPR

Then install it:

.. tab-set::
    :sync-group: tool

    .. tab-item:: uv (project)
        :sync: uv

        ``uv sync`` creates a ``.venv`` in the repository and installs PyEPR in
        editable mode, using the exact versions pinned in ``uv.lock``:

        .. code-block:: bash

            uv sync

        To also install the development (testing) and documentation tools:

        .. code-block:: bash

            uv sync --extra dev --extra doc

    .. tab-item:: uv (pip interface)
        :sync: uv-pip

        .. code-block:: bash

            uv venv
            source .venv/bin/activate      # Windows: .venv\Scripts\activate
            uv pip install -e ".[dev,doc]"

        Leave out ``[dev,doc]`` if you do not need the testing and
        documentation tools.

    .. tab-item:: pip
        :sync: pip

        .. code-block:: bash

            python -m venv .venv
            source .venv/bin/activate      # Windows: .venv\Scripts\activate
            python -m pip install -e ".[dev,doc]"

        Leave out ``[dev,doc]`` if you do not need the testing and
        documentation tools.

The ``-e`` (editable) flag means changes to the source code take effect without
reinstalling.

Verifying the installation
++++++++++++++++++++++++++

Check that PyEPR imports correctly and print its version:

.. code-block:: bash

    python -c "import pyepr; print(pyepr.__version__)"

or, inside a uv project:

.. code-block:: bash

    uv run python -c "import pyepr; print(pyepr.__version__)"

Requirements
++++++++++++

The following are installed automatically with PyEPR:

==========================================================  ================
Package                                                     Version
==========================================================  ================
Python                                                      >= 3.11
`NumPy <https://numpy.org>`_                                >= 2.1, < 3.0
`SciPy <https://scipy.org>`_                                >= 1.14.1
`Matplotlib <https://matplotlib.org>`_                      >= 3.9.2
`xarray <https://xarray.dev>`_                              >= 2025
`h5py <https://www.h5py.org>`_                              > 3.15
`h5netcdf <https://h5netcdf.org>`_                          > 1.4.0
`DeerLab <https://jeschkelab.github.io/DeerLab/>`_          >= 1.1.4, < 2.0
PyYAML                                                      >= 6.0.2, < 7.0
toml                                                        >= 0.10.2
psutil                                                      >= 7.1.3
requests                                                    >= 2.32
==========================================================  ================

Optional extras:

- ``dev``: testing tools (``pytest``, ``gitpython``).
- ``doc``: tools to build this documentation (``sphinx``, ``furo``,
  ``sphinx-gallery``, and related extensions).

Install them with ``uv sync --extra dev`` (or ``--extra doc``) from a source
checkout, or with ``pip install "pyepr-esr[dev]"``.

Hardware-specific requirements
++++++++++++++++++++++++++++++

You need these only to control a spectrometer. They cannot be installed from
PyPI and are not installed with PyEPR.

- **Bruker spectrometers (Xepr):** Xepr 2.9 or newer must be installed on the
  spectrometer computer, with its Python API enabled
  (*Processing → XeprAPI → Enable XeprAPI* in Xepr). PyEPR finds the
  ``XeprAPI`` module with ``Xepr --apipath``.
- **ETH AWG spectrometers:** requires MATLAB and the
  `MATLAB Engine API for Python <https://pypi.org/project/matlabengine/>`_,
  whose version must match your MATLAB release.
