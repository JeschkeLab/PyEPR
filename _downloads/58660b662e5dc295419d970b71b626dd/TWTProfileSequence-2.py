import matplotlib.pyplot as plt
import pyepr as epr
epr.set_waveform_precision(1)

seq = epr.TWTProfileSequence(
    B=12000, freq=34, reptime=3e3, averages=1, shots=10, dtp=2)
[tp, scale], data = seq.simulate(nu1=75, sat=0.5)
plt.pcolormesh(scale, tp, data.real, shading='auto')
plt.xlabel("Amplitude Scale")
plt.ylabel("Test Pulse Length (ns)")