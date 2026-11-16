# This file contains the object definitions for DEER Sequences
# 

import pyepr.pulses as p
from pyepr.sequences import Sequence
from warnings import warn
import numpy as np
import deerlab as dl
from pyepr.utils import val_in_us, add_phaseshift
from pyepr.classes import Parameter

class FourPulseDEERSequence(Sequence):
    """
    Class for the 4-pulse DEER sequence.
    """

    def __init__(self, *, B, freq,reptime,averages, shots,
                 tau1, tau2, step=8, start=None, dim=None,
                 phase_cycle='16-step',
                 **kwargs):
        r"""
        Initialize the 4-pulse DEER sequence.

        Sequence:
        [\pi/2 - \tau_1 - \pi - t \pi_{pump} - (\tau_2+\tau_1-t) - \pi - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau1 : int or float
            The first interpulse delay in us
        tau2 : int or float
            The second interpulse delay in us
        step : int or float
            The step size for the DEER time in ns, by default it is 8 ns
        start : int or float
            The starting DEER time in us, also known as tmin. If it is None then it is 180ns before tau2. By default it is None.
        dim : int
            The number of points in the DEER time, by default it is calculated from the step and start values, and stops at tau2
        phase_cycle : str|list, optional
            The phase cycling scheme to use. Either a pre-defined option 
            ('16-step' or 'DC'), can be supplied or a list of the phase cycles 
            for each pulse. Default is '16-step'. 

        Examples
        --------

        .. plot::
            :include-source:
        
            import pyepr as epr
            import numpy as np
            from pyepr.sequences.DEER import FourPulseDEERSequence
            N_gyro = 0.002803632236095
            freq = 34

            exc_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi/2)
            ref_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi)
            pump_pulse = epr.ChirpPulse(tp=200,init_freq=-0.10,final_freq=-0.05,flipangle=np.pi)
            seq = FourPulseDEERSequence(
                    B=freq/N_gyro,
                    freq=freq,
                    reptime=3e3,
                    averages=1,
                    shots=10,
                    tau1=0.4,
                    tau2=4,
                    step=8,
                    start=0.28,
                    phase_cycle='16-step',
                    exc_pulse=exc_pulse,
                    ref_pulse=ref_pulse,
                    pump_pulse=pump_pulse
                )

            seq.draw()
            print(seq)
        """
        name = kwargs.pop('name',"FourPulseDEER")

        self._base_pulses = {}
        self._base_pulses['exc'] = kwargs.pop('exc_pulse', p.RectPulse(tp=16, freq=0, flipangle=np.pi/2))
        self._base_pulses['pump'] = kwargs.pop('pump_pulse', p.RectPulse(tp=32, freq=-0.1, flipangle=np.pi))
        self._base_pulses['ref'] = kwargs.pop('ref_pulse', p.RectPulse(tp=32, freq=0, flipangle=np.pi))
        self._base_pulses['det'] = kwargs.pop('det_pulse', p.Detection(tp=256,freq=self._base_pulses['exc'].freq.value))

        if self._base_pulses['det'].freq.value != self._base_pulses['exc'].freq.value:
            warn("The detection pulse frequency is not equal to the excitation pulse frequency. This may lead to unexpected results.")


        self.tau1: Parameter = Parameter(name="tau1", value=tau1*1e3, unit="ns", description="Interpulse Delay 1", virtual=True)
        """The first interpulse delay parameter in ns."""
        self.tau2: Parameter = Parameter(name="tau2", value=tau2*1e3, unit="ns", description="Interpulse Delay 2", virtual=True)
        """The second interpulse delay parameter in ns."""
        if start is None:
            start = tau1*1e3 - 180 # Start 180ns before tau2
        else:
            start *= 1e3 # Convert to ns
        if dim is None:
            dim = int(tau2*1e3 / step) + 1

        self.t: Parameter = Parameter(name="t", value=start, step=step, dim=dim, unit="ns", description="DEER Time", virtual=True)
        """The pump pulse position parameter in ns."""
        super().__init__(name=name, B=B, freq=freq, reptime=reptime, averages=averages, shots=shots)

        self._build_sequence()
        self._build_pcyc(phase_cycle)

        self.evolution([self.t])

    def _build_sequence(self):

        self.addPulse(self._base_pulses['exc'].copy(t=0, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['ref'].copy(t=self.tau1, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['pump'].copy(t=self.tau1+self.t, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['ref'].copy(t=2*self.tau1+self.tau2, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['det'].copy(t=2*(self.tau1+self.tau2)+self._base_pulses['ref'].tp.value, pcyc={'phases':[0], 'dets':[1]}))


    def _build_pcyc(self, phase_cycle):

        if isinstance(phase_cycle,list):
            if len(phase_cycle) != 4:
                raise ValueError("Phase cycle list must have 4 elements, one for each pulse.")
            self._pcyc = phase_cycle
            self.pcyc_name = "Custom Phase Cycle"

        elif phase_cycle == '16-step':
            self._pcyc = [
                [[0],[1]], # No phase cycling for the excitation pulse
                [[0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1]],
                [[0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1]],
                [[0],[1]],
            ]
            self.pcyc_name = "16-step"
        elif phase_cycle == '2-step': # Used only when the pump pulse is 
            self._pcyc = [
                [[0,np.pi],[1,-1]], # No phase cycling for the excitation pulse
                [[0],[1]],
                [[0],[1]],
                [[0],[1]],
            ]
            self.pcyc_name = "2-step"
        else:
            raise ValueError("Invalid phase cycle option. Choose from '16-step', '2-step', or provide a list of 4 phase cycles.")

        for i in range(len(self.pulses)-1):
            self.pulses[i]._addPhaseCycle(*self._pcyc[i])

        self._buildPhaseCycle()


    def simulate(self,**kwargs):
        tau1 = val_in_us(self.tau1)
        tau2 = val_in_us(self.tau2)
        t = val_in_us(self.t)

        experimentInfo = dl.ex_4pdeer(tau1=tau1,tau2=tau2,pathways=[1,2,3])
        reftimes = dict(zip(["reftime1","reftime2","reftime3"],experimentInfo.reftimes(tau1,tau2)))
        mod_depths = {"lam1":0.4, "lam2":0.1, "lam3":0.2}

        r = np.linspace(0.5,10,100)
        rmean = kwargs.pop('rmean', 4.0)
        rstd = kwargs.pop('rstd', 0.75)
        conc = kwargs.pop('conc', 100)

        Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experimentInfo)
        Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)

        Vsim = add_phaseshift(Vsim, 0.05)
        return t, Vsim

    
class FivePulseDEERSequence(Sequence):
                 
    pass


class ThreePulseDEERSequence(Sequence):

    pass