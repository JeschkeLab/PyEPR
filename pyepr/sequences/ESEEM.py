import pyepr.pulses as p
from pyepr.sequences import Sequence
from pyepr.classes import Parameter
from warnings import warn
import numpy as np


class ThreePulseESEEMSequence(Sequence):

    def __init__(self, *, B, freq, reptime, averages, shots, tau=300, start=500, step=40, dim=200, **kwargs) -> None:
        """
        Creates a Three Pulse ESEEM sequence object. 
        In the form of [\pi/2 - \tau - \pi/2 - T - \pi/2 - echo]

        Parameters
        ----------
        B : float
            The B0 field in Gauss.
        freq : float
            The frequency in GHz.
        reptime : float
            The shot repetition time in us.
        averages : int
            The number of scans.
        shots : int
            The number of shots per point.
        tau : float, optional
            The first interpulse delay in us, by default 300.
        start : float, optional
            The starting mixing time in ns, by default 500.
        step : float, optional
            The step size for the mixing time in ns, by default 40.
        dim : int, optional
            The number of points in the mixing time, by default 200.
        pi2_pulse : Pulse, optional
            The pi/2 pulse to use, by default None. If None, a RectPulse will be used.
        pi_pulse : Pulse, optional
            The pi pulse to use, by default None. If None, a RectPulse will be used.


        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)
            seq = epr.ThreePulseESEEMSequence(
                B=12000, freq=34, reptime=1000, averages=1, shots=10,
                tau=300, start=500, step=40, dim=200)
            seq.draw() 
        
        """

        self.tau: Parameter = Parameter(name="tau", value=start,step=step,dim=dim, unit="ns", description="1st Inter Pulse Delay",virtual=True)
        """First Interpulse Delay parameter in ns."""
        
        self.T: Parameter = Parameter(
            name='T',
            description = 'Mixing Time',
            unit = 'ns',
            value=start,
            step = step,
            dim = dim,
            virtual=True)
        """1st Mixing Time parameter in ns."""

        super().__init__(name="ThreePulseESEEMSequence",
                         B=B, freq=freq, reptime=reptime, averages=averages, shots=shots, **kwargs)

        self._build_sequence()

    def _build_sequence(self):
        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=0, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(p.RectPulse(  # pi/2
                t=0, tp=16, freq=0, flipangle=np.pi/2, 
                pcyc={"phases":[0, np.pi], "dets": [1, -1]}
            ))

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=self.tau))
        else:
            self.addPulse(p.RectPulse(  # pi/2
            t=self.tau, tp=32, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=self.tau+self.T))
        else:
            self.addPulse(p.RectPulse(  # pi/2
            t=self.tau+self.T, tp=32, freq=0, flipangle=np.pi
            ))

        self.addPulse(p.Detection(t=self.tau+self.T+self.tau, tp=512))

        self.evolution([self.T])

    def _simulate(self):
        raise NotImplementedError("Simulation for Three Pulse ESEEM is not implemented yet.")

    
class HYSCORESequence(Sequence):
    def __init__(self, *, B, freq, reptime, averages, shots, tau=300, start=300, step=16, dim=128, **kwargs) -> None:
        """
        Creates a HYSCORE sequence object.
        In the form of [\pi/2 - \tau - \pi/2 - T1 - \pi - T2 - \pi/2 - echo]

        Parameters
        ----------
        B : float
            The B0 field in Gauss.
        freq : float
            The frequency in GHz.
        reptime : float
            The shot repetition time in us.
        averages : int
            The number of scans.
        shots : int
            The number of shots per point.
        tau : float, optional
            The first interpulse delay in us, by default 300.
        start : float, optional
            The starting mixing time in ns, by default 300.
        step : float, optional
            The step size for the mixing time in ns, by default 16.
        dim : int, optional
            The number of points in the mixing time, by default 128.

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.HYSCORESequence(
                B=12000, freq=34, reptime=1000, averages=1, shots=10,
                tau=300, start=300, step=16, dim=128)
            seq.draw()

        """
        
        self.tau: Parameter = Parameter(
            name="tau", value=start,step=step, dim=dim, unit="ns",
            description="1st Inter Pulse Delay",virtual=True)
        """First Interpulse Delay parameter in ns."""

        self.T1: Parameter = Parameter(
            name='T1',
            description = '1st Mixing Time',
            unit = 'ns',
            value=start,
            step = step,
            dim = dim,
            virtual=True)
        """1st Mixing Time parameter in ns."""
        
        self.T2: Parameter = Parameter(
            name='T2',
            description = '2nd Mixing Time',
            unit = 'ns',
            value=start,
            step = step,
            dim = dim,
            virtual=True)
        """2nd Mixing Time parameter in ns."""

        super().__init__(
            name = "HYSCORESequence",
            B=B, freq=freq, reptime=reptime, averages=averages, shots=shots, **kwargs)

        self._build_sequence()

    def _build_sequence(self):

        if hasattr(self, "pi2_pulse"):
            p90 = self.pi2_pulse.copy(t=0, pcyc={'phases':[0], 'dets':[1]})
        else:
            p90 = p.RectPulse(  # pi/2
                t=0, tp=16, freq=0, flipangle=np.pi/2, 
                pcyc={"phases":[0], "dets": [1]}
            )

        if hasattr(self, "pi_pulse"):
            p180 = self.pi_pulse.copy(t=0, pcyc={'phases':[0], 'dets':[1]})
        else:
            p180 = p.RectPulse(  # pi
                t=0, tp=32, freq=0, flipangle=np.pi, 
                pcyc={"phases":[0], "dets": [1]}
            )

        det_event = p.Detection(t=0, tp=512)


        self.addPulse(p90.copy(t=0, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(p90.copy(t=self.tau, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(p180.copy(t=self.tau+self.T1,pcyc={'phases':[0,np.pi/2,np.pi,3*np.pi/2,], 'dets':[1,-1,1,-1]}))
        self.addPulse(p90.copy(t=self.tau+self.T1+self.T2,pcyc={'phases':[0], 'dets':[1]}))

        self.addPulse(det_event.copy(t=2*self.tau+self.T1+self.T2))
        self._buildPhaseCycle()
        self.evolution([self.T1,self.T2])

    def simulate(self):
        raise NotImplementedError("Simulation for HYSCORE is not implemented yet.")