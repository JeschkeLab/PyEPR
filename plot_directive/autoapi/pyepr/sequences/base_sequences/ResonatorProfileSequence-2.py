import matplotlib.pyplot as plt
import pyepr as epr
import numpy as np
epr.set_waveform_precision(1)

N_gyro = 0.002803
freq = 34
seq = epr.ResonatorProfileSequence(
    B=freq/N_gyro, freq=freq, reptime=3e3, averages=1, shots=10,
    fwidth=0.3, fstep=0.02, dtp=1, dim=100)
[tp, f], data = seq.simulate(Q=100, nu1=75)
plt.pcolormesh(f, tp, data.real, shading='auto')
plt.xlabel("Frequency (GHz)")
plt.ylabel("Test Pulse Length (ns)")