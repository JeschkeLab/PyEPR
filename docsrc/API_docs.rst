.. rst-class:: api-docs

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
   pyepr.AmplifierLinearityAnalysis

Sequences
~~~~~~~~~

Setup and Calibration Sequences
+++++++++++++++++++++++++++++++
.. _SetupSequences:
.. autoapisummary::

   pyepr.sequences.HahnEchoSequence
   pyepr.sequences.FieldSweepSequence   
   pyepr.sequences.ResonatorProfileSequence
   pyepr.sequences.TWTProfileSequence

Relxation Sequences
+++++++++++++++++++++
.. _RelaxationSequences:
.. autoapisummary::

    pyepr.sequences.HahnEchoRelaxationSequence
    pyepr.sequences.CarrPurcellSequence
    pyepr.sequences.T1InversionRecoverySequence
    pyepr.sequences.ReptimeScanSequence

Dipolar EPR Sequences
+++++++++++++++++++++
.. _DipolarSequences:

.. autoapisummary::

   pyepr.sequences.DEER.ThreePulseDEERSequence
   pyepr.sequences.DEER.FourPulseDEERSequence
   pyepr.sequences.DEER.FivePulseDEERSequence

Hyperfine Spectroscopy Sequences
++++++++++++++++++++++++++++++++
.. _HyperfineSequences:

.. autoapisummary::

   pyepr.sequences.ESEEM.ThreePulseESEEMSequence
   pyepr.sequences.ESEEM.HYSCORESequence

Pulses
~~~~~~
.. _Pulses:

.. autoapisummary::

   pyepr.pulses.Pulse
   pyepr.pulses.Detection
   pyepr.pulses.RectPulse
   pyepr.pulses.GaussianPulse
   pyepr.pulses.ChirpPulse
   pyepr.pulses.HSPulse
   pyepr.pulses.SincPulse
   pyepr.pulses.WURSTPulse

Termination Criteria
~~~~~~~~~~~~~~~~~~~~

Termination criteria are used for stopping sequences under specific conditions,
such as when a certain signal-to-noise ratio is reached or when a specific time has elapsed. 
They are used through the :py:meth:`interface.terminate_at <pyepr.hardware.BaseInterfaceClass.Interface.terminate_at>` method.

.. autoapisummary::

    pyepr.criteria.Criteria
    pyepr.criteria.TimeCriteria
    pyepr.criteria.SNRCriteria


Interfaces
~~~~~~~~~~
.. _Interfaces:

.. autoapisummary::

    pyepr.hardware.BaseInterfaceClass.Interface
    pyepr.hardware.Bruker_AWG.BrukerAWG
    pyepr.hardware.Bruker_MPFU.BrukerMPFU
    pyepr.hardware.ETH_awg.ETH_awg_interface
    pyepr.hardware.PyEPR_control.PyEPRControlInterface
    pyepr.hardware.XeprAPI_link.XeprAPILink

Utilities
~~~~~~~~~

.. autoapisummary::

    pyepr.dataset.EPRAccessor

Functions
----------------

I/O
~~~

.. autoapisummary::
    :nosignatures:

    pyepr.tools.eprload
    pyepr.utils.save_file
    pyepr.dataset.create_dataset_from_sequence
    pyepr.dataset.create_dataset_from_axes
    pyepr.dataset.create_dataset_from_bruker
    pyepr.dataset.downconvert_dataset

Utilities
~~~~~~~~~

.. autoapisummary::
    :nosignatures:

    pyepr.utils.transpose_dict_of_list
    pyepr.utils.transpose_list_of_dicts
    pyepr.utils.round_step
    pyepr.utils.gcd
    pyepr.utils.sop
    pyepr.hardware.Bruker_tools.write_pulsespel_file