Tutorial
========

It is recommended that most users should use PyEPR inside a Jupyter notebook, unless PyEPR is being used as a backend for another program.

If you have never used Jupyter notebooks before, we recommend that you first familiarize yourself with the basics of Jupyter notebooks.
There are many tutorials available online, such as the one provided by the Jupyter project. <https://docs.jupyter.org/en/latest/>

To get started, we will first import the PyEPR package.

.. code-block:: python

    import pyepr as epr

This will import the PyEPR package, except for the hardware modules. Importing of hardware modules is included hardware-control tutorial.

Tutorial Overview
-----------------
The following tutorials are available:

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: Loading and Analysing Data
        :link: tutorial_loading
        :link-type: doc

        Learn how to load data from the PyEPR format and other common formats, and perform basic analysis with PyEPR.

    .. grid-item-card:: Sequencer
        :link: tutorial_sequencer
        :link-type: doc

        Explore how to create pulse sequences, and use default sequences provided by PyEPR.

    .. grid-item-card:: Hardware Control
        :link: tutorial_hardware
        :link-type: doc

        Interface with physical hardware devices.

    .. grid-item-card:: Examples
        :link: auto_examples/index
        :link-type: doc

        View practical examples and use cases.


.. toctree::
    :maxdepth: 1
    :hidden:
    
    tutorial_loading.md
    tutorial_sequencer.md
    tutorial_hardware.md
    auto_examples/index.rst