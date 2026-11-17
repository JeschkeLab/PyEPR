import matplotlib.pyplot as plt
import pyepr as epr
epr.set_waveform_precision(1)

seq = epr.T1InversionRecoverySequence(
    B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=500, step=500, dim=100)
x, data = seq.simulate(T1=1e4)
plt.plot(x, data.real)
plt.xlabel("Inversion Recovery Delay (ns)")
plt.ylabel("Signal Intensity (a.u.)")