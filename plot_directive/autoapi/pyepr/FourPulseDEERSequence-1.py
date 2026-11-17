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