import pyepr as epr
import numpy as np
epr.set_waveform_precision(1)

N_gyro = 0.002803
freq = 34
pi2_pulse = epr.RectPulse(tp=32,freq=0,flipangle=np.pi/2)
pi_pulse = epr.RectPulse(tp=32,freq=0,flipangle=np.pi)

seq = epr.ResonatorProfileSequence(
    B=freq/N_gyro,
    freq=freq,
    reptime=3e3,
    averages=1,
    shots=10,
    pi2_pulse=pi2_pulse,
    pi_pulse=pi_pulse,
    fwidth=0.3,
    fstep=0.02,
    dtp=1,
    dim=40)

seq.draw()