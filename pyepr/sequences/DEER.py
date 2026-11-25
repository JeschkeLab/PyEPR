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
        reptime : float
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
        phase_cycle : str or list, optional
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
        r"""
        Simulate the 4-pulse DEER sequence using DeerLab. 

        Optional Parameters
        -------------------
        rmean : float
            The mean distance of the spin pair distribution in nm. Default is 4.0 nm
        rstd : float
            The standard deviation of the spin pair distribution in nm. Default is 0.75 nm
        conc : float
            The concentration of the spin pair in uM. Default is 100 uM
        
        Returns
        -------
        t : np.ndarray
            The DEER time axis in us
        Vsim : np.ndarray
            The simulated DEER signal
        
        Examples
        --------

        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np
            import matplotlib.pyplot as plt
            N_gyro = 0.002803632236095
            freq = 34

            exc_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi/2)
            ref_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi)
            pump_pulse = epr.ChirpPulse(tp=200,init_freq=-0.10,final_freq=-0.05,flipangle=np.pi)
            seq = epr.FourPulseDEERSequence(
                    B=freq/N_gyro,
                    freq=freq,
                    reptime=3e3,
                    averages=1,
                    shots=10,
                    tau1=0.4,
                    tau2=3,
                    step=8,
                    start=0.1,
                    phase_cycle='16-step',
                    exc_pulse=exc_pulse,
                    ref_pulse=ref_pulse,
                    pump_pulse=pump_pulse
                )

            t, Vsim = seq.simulate(rmean=4.0, rstd=0.75, conc=100)
            plt.plot(t, Vsim)
            plt.xlabel(r'Time ($\mu$s)')
            plt.ylabel(r'Signal (a.u.)')
        """

        tau1 = val_in_us(self.tau1)
        tau2 = val_in_us(self.tau2)
        t = val_in_us(self.t)

        experimentInfo = dl.ex_4pdeer(tau1=tau1,tau2=tau2,pathways=[1,2,3])
        reftimes = dict(zip(["reftime1","reftime2","reftime3"],experimentInfo.reftimes(tau1,tau2)))
        mod_depths = {"lam1":0.4, "lam2":0.05, "lam3":0.05}

        r = np.linspace(0.5,10,100)
        rmean = kwargs.pop('rmean', 4.0)
        rstd = kwargs.pop('rstd', 0.75)
        conc = kwargs.pop('conc', 100)

        Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experimentInfo)
        Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)

        Vsim = add_phaseshift(Vsim, 0.05)
        return t, Vsim

    
class FivePulseDEERSequence(Sequence):

    def __init__(self, *, B, freq,reptime,averages, shots,
        tau1, tau2, tau3, step=8, start=None, dim=None,
        phase_cycle='16-step',
        **kwargs):
        r"""
        Initialize the 5-pulse DEER sequence.

        Sequence:
        [\pi/2 - \tau_1-\tau_3 \pi_{pump} - \tau3 \pi - t \pi_{pump} - (\tau_2+\tau_1-t) - \pi - \tau2 - echo]   

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : float
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau1 : int or float
            The first interpulse delay in us
        tau2 : int or float
            The second interpulse delay in us
        tau3 : int or float
            The third interpulse delay in us
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
            from pyepr.sequences.DEER import FivePulseDEERSequence
            N_gyro = 0.002803632236095
            freq = 34

            exc_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi/2)
            ref_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi)
            pump_pulse = epr.ChirpPulse(tp=200,init_freq=-0.10,final_freq=-0.05,flipangle=np.pi)
            seq = FivePulseDEERSequence(
                    B=freq/N_gyro,
                    freq=freq,
                    reptime=3e3,
                    averages=1,
                    shots=10,
                    tau1=3,
                    tau2=3,
                    tau3=0.3,
                    step=8,
                    start=0.1,
                    phase_cycle='16-step',
                    exc_pulse=exc_pulse,
                    ref_pulse=ref_pulse,
                    pump_pulse=pump_pulse
                )

            seq.draw()
            print(seq)

        """

        name = kwargs.pop('name',"FivePulseDEER")

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
        self.tau3: Parameter = Parameter(name="tau3", value=tau3*1e3, unit="ns", description="Interpulse Delay 3", virtual=True)
        """The third interpulse delay parameter in ns."""
        if start is None:
            start = tau3*1e3 - 180 # Start 180ns before tau2
        else:
            start *= 1e3 # Convert to ns
        if dim is None:
            dim = int((tau1+tau2)*1e3 / step) + 1

        self.t: Parameter = Parameter(name="t", value=start, step=step, dim=dim, unit="ns", description="DEER Time", virtual=True)
        """The pump pulse position parameter in ns."""
        super().__init__(name=name, B=B, freq=freq, reptime=reptime, averages=averages, shots=shots)

        self._build_sequence()
        self._build_pcyc(phase_cycle)

        self.evolution([self.t])

    def _build_sequence(self):

        self.addPulse(self._base_pulses['exc'].copy(t=0, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['pump'].copy(t=self.tau1-self.tau3, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['ref'].copy(t=self.tau1, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['pump'].copy(t=self.tau1+self.t, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['ref'].copy(t=2*self.tau1+self.tau2, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['det'].copy(t=2*(self.tau1+self.tau2), pcyc={'phases':[0], 'dets':[1]}))

    def _build_pcyc(self, phase_cycle):
        
        if isinstance(phase_cycle,list):
            if len(phase_cycle) != 5:
                raise ValueError("Phase cycle list must have 5 elements, one for each pulse.")
            self._pcyc = phase_cycle
            self.pcyc_name = "Custom Phase Cycle"

        elif phase_cycle == '16-step':
            self._pcyc = [
                [[0],[1]], # No phase cycling for the excitation pulse
                [[0],[1]],
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
                [[0],[1]],
            ]
            self.pcyc_name = "2-step"
        else:
            raise ValueError("Invalid phase cycle option. Choose from '16-step', '2-step', or provide a list of 5 phase cycles.")

        for i in range(len(self.pulses)-1):
            self.pulses[i]._addPhaseCycle(*self._pcyc[i])

        self._buildPhaseCycle()

    def simulate(self,**kwargs):
        r"""
        Simulate the 5-pulse DEER sequence using DeerLab. 

        Optional Parameters
        -------------------
        rmean : float
            The mean distance of the spin pair distribution in nm. Default is 4.0 nm
        rstd : float
            The standard deviation of the spin pair distribution in nm. Default is 0.75 nm
        conc : float
            The concentration of the spin pair in uM. Default is 100 uM
        
        Returns
        -------
        t : np.ndarray
            The DEER time axis in us
        Vsim : np.ndarray
            The simulated DEER signal
        
        Examples
        --------

        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np
            import matplotlib.pyplot as plt
            from pyepr.sequences.DEER import FivePulseDEERSequence
            N_gyro = 0.002803632236095
            freq = 34

            exc_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi/2)
            ref_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi)
            pump_pulse = epr.ChirpPulse(tp=200,init_freq=-0.10,final_freq=-0.05,flipangle=np.pi)
            seq = FivePulseDEERSequence(
                    B=freq/N_gyro,
                    freq=freq,
                    reptime=3e3,
                    averages=1,
                    shots=10,
                    tau1=3,
                    tau2=3,
                    tau3=0.3,
                    step=8,
                    start=0.1,
                    phase_cycle='16-step',
                    exc_pulse=exc_pulse,
                    ref_pulse=ref_pulse,
                    pump_pulse=pump_pulse
                )

            t, Vsim = seq.simulate(rmean=4.0, rstd=0.75, conc=100)
            plt.plot(t, Vsim)
            plt.xlabel(r'Time ($\mu$s)')
            plt.ylabel(r'Signal (a.u.)')
        """

        tau1 = val_in_us(self.tau1)
        tau2 = val_in_us(self.tau2)
        tau3 = val_in_us(self.tau3)
        t = val_in_us(self.t)

        experimentInfo = dl.ex_fwd5pdeer(tau1=tau1,tau2=tau2,tau3=tau3,pathways=[1,5])
        reftimes = dict(zip(["reftime1","reftime5"],experimentInfo.reftimes(tau1,tau2,tau3)))
        mod_depths = {"lam1":0.4, "lam5":0.1}

        r = np.linspace(0.5,10,100)
        rmean = kwargs.pop('rmean', 4.0)
        rstd = kwargs.pop('rstd', 0.75)
        conc = kwargs.pop('conc', 100)

        Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experimentInfo)
        Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)

        Vsim = add_phaseshift(Vsim, 0.05)
        return t, Vsim





class ThreePulseDEERSequence(Sequence):

    def __init__(self, *, B, freq, reptime, averages, shots,
                 tau1, step=8, start=None, dim=None,
                 phase_cycle='16-step',
                 **kwargs):
        r"""
        Initialize the 3-pulse DEER sequence.

        Sequence:
        [\pi/2 - t - \pi_{pump} - (\tau_1-t) - \pi - \tau_1 - echo]

        The DEER time `t` is the delay between the excitation pulse and the
        pump pulse. The pump pulse is swept across the first interpulse delay,
        the main dipolar pathway therefore refocuses at `t=0`.

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
            The interpulse delay in us
        step : int or float
            The step size for the DEER time in ns, by default it is 8 ns
        start : int or float
            The starting DEER time in us. If it is None then the pump pulse
            starts immediately after the excitation pulse, so that the two
            pulses do not overlap. By default it is None.
        dim : int
            The number of points in the DEER time, by default it is calculated
            from the step and start values, and stops before the pump pulse
            overlaps with the refocusing pulse.
        phase_cycle : str|list, optional
            The phase cycling scheme to use. Either a pre-defined option
            ('8-step' or '2-step'), can be supplied or a list of the phase
            cycles for each pulse. Default is '8-step'.

        Examples
        --------

        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np
            N_gyro = 0.002803632236095
            freq = 34

            exc_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi/2)
            ref_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi)
            pump_pulse = epr.ChirpPulse(tp=200,init_freq=-0.10,final_freq=-0.05,flipangle=np.pi)
            seq = epr.ThreePulseDEERSequence(
                    B=freq/N_gyro,
                    freq=freq,
                    reptime=3e3,
                    averages=1,
                    shots=10,
                    tau1=3,
                    step=8,
                    phase_cycle='8-step',
                    exc_pulse=exc_pulse,
                    ref_pulse=ref_pulse,
                    pump_pulse=pump_pulse
                )

            seq.draw()
            print(seq)
        """
        name = kwargs.pop('name', "ThreePulseDEER")

        self._base_pulses = {}
        self._base_pulses['exc'] = kwargs.pop('exc_pulse', p.RectPulse(tp=16, freq=0, flipangle=np.pi/2))
        self._base_pulses['pump'] = kwargs.pop('pump_pulse', p.RectPulse(tp=32, freq=-0.1, flipangle=np.pi))
        self._base_pulses['ref'] = kwargs.pop('ref_pulse', p.RectPulse(tp=32, freq=0, flipangle=np.pi))
        self._base_pulses['det'] = kwargs.pop('det_pulse', p.Detection(tp=256, freq=self._base_pulses['exc'].freq.value))

        if self._base_pulses['det'].freq.value != self._base_pulses['exc'].freq.value:
            warn("The detection pulse frequency is not equal to the excitation pulse frequency. This may lead to unexpected results.")

        self.tau1: Parameter = Parameter(name="tau1", value=tau1*1e3, unit="ns", description="Interpulse Delay 1", virtual=True)
        """The interpulse delay parameter in ns."""

        exc_tp = self._base_pulses['exc'].tp.value
        pump_tp = self._base_pulses['pump'].tp.value
        ref_tp = self._base_pulses['ref'].tp.value

        if start is None:
            start = exc_tp+5  # Start 5ns after the excitation pulse
        else:
            start *= 1e3  # Convert to ns
        if dim is None:
            stop = tau1*1e3 - (pump_tp + ref_tp) / 2  # Pump directly before the refocusing pulse
            dim = int((stop - start) / step) + 1
        if dim < 1:
            raise ValueError("tau1 is too short for the given pulses, start and step values.")

        self.t: Parameter = Parameter(name="t", value=start, step=step, dim=dim, unit="ns", description="DEER Time", virtual=True)
        """The pump pulse position parameter in ns."""
        super().__init__(name=name, B=B, freq=freq, reptime=reptime, averages=averages, shots=shots)

        self._build_sequence()
        self._build_pcyc(phase_cycle)

        self.evolution([self.t])

    def _build_sequence(self):

        self.addPulse(self._base_pulses['exc'].copy(t=0, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['pump'].copy(t=self.t, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['ref'].copy(t=self.tau1, pcyc={'phases':[0], 'dets':[1]}))
        self.addPulse(self._base_pulses['det'].copy(t=2*self.tau1, pcyc={'phases':[0], 'dets':[1]}))

    def _build_pcyc(self, phase_cycle):

        if isinstance(phase_cycle, list):
            if len(phase_cycle) != 3:
                raise ValueError("Phase cycle list must have 3 elements, one for each pulse.")
            self._pcyc = phase_cycle
            self.pcyc_name = "Custom Phase Cycle"

        elif phase_cycle == '8-step':
            self._pcyc = [
                [[0, np.pi],[1,-1]], # No phase cycling for the excitation pulse
                [[0],[1]],
                [[0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1]],
            ]
            self.pcyc_name = "8-step"
        elif phase_cycle == '2-step':
            self._pcyc = [
                [[0, np.pi], [1, -1]],
                [[0],[1]],
                [[0],[1]],
            ]
            self.pcyc_name = "2-step"
        else:
            raise ValueError("Invalid phase cycle option. Choose from '16-step', '2-step', or provide a list of 3 phase cycles.")

        for i in range(len(self.pulses)-1):
            self.pulses[i]._addPhaseCycle(*self._pcyc[i])

        self._buildPhaseCycle()

    def simulate(self, **kwargs):
        r"""
        Simulate the 3-pulse DEER sequence using DeerLab.

        Optional Parameters
        -------------------
        rmean : float
            The mean distance of the spin pair distribution in nm. Default is 4.0 nm
        rstd : float
            The standard deviation of the spin pair distribution in nm. Default is 0.75 nm
        conc : float
            The concentration of the spin pair in uM. Default is 100 uM

        Returns
        -------
        t : np.ndarray
            The DEER time axis in us
        Vsim : np.ndarray
            The simulated DEER signal

        Examples
        --------

        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np
            import matplotlib.pyplot as plt
            from pyepr.sequences.DEER import ThreePulseDEERSequence
            N_gyro = 0.002803632236095
            freq = 34

            exc_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi/2)
            ref_pulse = epr.GaussianPulse(tp=32,freq=0,flipangle=np.pi)
            pump_pulse = epr.ChirpPulse(tp=200,init_freq=-0.10,final_freq=-0.05,flipangle=np.pi)
            seq = ThreePulseDEERSequence(
                    B=freq/N_gyro,
                    freq=freq,
                    reptime=3e3,
                    averages=1,
                    shots=10,
                    tau1=3,
                    step=8,
                    phase_cycle='8-step',
                    exc_pulse=exc_pulse,
                    ref_pulse=ref_pulse,
                    pump_pulse=pump_pulse
                )

            t, Vsim = seq.simulate(rmean=4.0, rstd=0.75, conc=100)
            plt.plot(t, Vsim)
            plt.xlabel(r'Time ($\\mu$s)')
            plt.ylabel(r'Signal (a.u.)')
        """

        tau1 = val_in_us(self.tau1)
        t = val_in_us(self.t)

        # DeerLab pathway 1 refocuses at t=tau, pathway 2 at t=0 (main pathway)
        experimentInfo = dl.ex_3pdeer(tau=tau1, pathways=[1,2])
        reftimes = dict(zip(["reftime1","reftime2"], experimentInfo.reftimes(tau1)))
        mod_depths = {"lam1":0.05, "lam2":0.4}

        r = np.linspace(0.5,10,100)
        rmean = kwargs.pop('rmean', 4.0)
        rstd = kwargs.pop('rstd', 0.75)
        conc = kwargs.pop('conc', 100)

        Vmodel = dl.dipolarmodel(t, r, Pmodel=dl.dd_gauss, experiment=experimentInfo)
        Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)

        Vsim = add_phaseshift(Vsim, 0.05)
        return t, Vsim