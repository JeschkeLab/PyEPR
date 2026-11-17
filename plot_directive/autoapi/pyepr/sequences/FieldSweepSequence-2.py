import matplotlib.pyplot as plt
import pyepr as epr
epr.set_waveform_precision(1)

N_gyro = 0.002803632236095
freq = 34
seq = epr.FieldSweepSequence(
    B=freq/N_gyro, Bwidth=200, freq=freq, reptime=1e3, averages=1, shots=10)
B, data = seq.simulate()
plt.plot(B, data.real)
plt.xlabel("Magnetic Field (G)")
plt.ylabel("Signal Intensity (a.u.)")