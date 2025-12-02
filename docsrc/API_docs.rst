API doc
=================
Classes
----------------

Main Classes
~~~~~~~~~~~~

.. autoapisummary::
    
    pyepr.pulses.Pulse
    pyepr.criteria.Criteria
    pyepr.classes.Parameter
    pyepr.sequences.Sequence
    
Analysis Modules
~~~~~~~~~~~~~~~~

.. autoapisummary::

   pyepr.FieldSweepAnalysis
   pyepr.HahnEchoRelaxationAnalysis
   pyepr.ResonatorProfileAnalysis
   pyepr.CarrPurcellAnalysis
   pyepr.ReptimeAnalysis

Sequences
~~~~~~~~~
.. _Sequences:
.. autoapisummary::

   pyepr.sequences.HahnEchoSequence
   pyepr.sequences.HahnEchoRelaxationSequence
   pyepr.sequences.FieldSweepSequence
   pyepr.sequences.ReptimeScan
   pyepr.sequences.CarrPurcellSequence
   pyepr.sequences.ResonatorProfileSequence
   pyepr.sequences.TWTProfileSequence

Pulses
~~~~~~
.. _Pulses:

.. autoapisummary::

   pyepr.pulses.Pulse
   pyepr.pulses.Detection
   pyepr.pulses.RectPulse
   pyepr.pulses.GaussianPulse
   pyepr.pulses.HSPulse
   pyepr.pulses.ChirpPulse
   pyepr.pulses.SincPulse

Termination Criteria
~~~~~~~~~~~~~~~~~~~~

.. autoapisummary::

    pyepr.criteria.Criteria
    pyepr.criteria.TimeCriteria
    pyepr.criteria.SNRCriteria

Utilities
~~~~~~~~~

.. autoapisummary::

    pyepr.dataset.EPRAccessor

Interfaces
~~~~~~~~~~
.. _Interfaces:

.. autoapisummary::
    pyepr.classes.Interface
    pyepr.hardware.Bruker_AWG.BrukerAWG
    pyepr.hardware.Bruker_MPFU.BrukerMPFU
    pyepr.hardware.XeprAPI_link.XeprAPILink
    pyepr.hardware.ETH_awg.ETH_awg_interface

Functions
----------------

I/O
~~~

.. autoapisummary::
    pyepr.tools.eprload
    pyepr.utils.save_file
    pyepr.dataset.create_dataset_from_sequence
    pyepr.dataset.create_dataset_from_axes
    pyepr.dataset.create_dataset_from_bruker

Utilities
~~~~~~~~~

.. autoapisummary::
    pyepr.utils.transpose_dict_of_list
    pyepr.utils.transpose_list_of_dicts
    pyepr.utils.round_step
    pyepr.utils.gcd
    pyepr.utils.sop
    pyepr.hardware.Bruker_tools.write_pulsespel_file