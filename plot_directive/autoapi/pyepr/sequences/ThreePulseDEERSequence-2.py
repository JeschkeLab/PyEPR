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