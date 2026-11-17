import matplotlib.pyplot as plt
import pyepr as epr
epr.set_waveform_precision(1)

seq = epr.CarrPurcellSequence(B=12000, freq=34, reptime=1000, averages=1, shots=10, n=4)
x,data = seq.simulate(T_CP=2.5e3, e=1.8)
plt.plot(x,data)
plt.xlabel("Interpulse Delay (ns)")
plt.ylabel("Signal Intensity (a.u.)")