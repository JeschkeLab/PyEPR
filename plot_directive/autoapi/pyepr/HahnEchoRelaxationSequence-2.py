import matplotlib.pyplot as plt
import pyepr as epr
epr.set_waveform_precision(1)

seq = epr.HahnEchoRelaxationSequence(
    B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=300, step=40, dim=200)
x, data = seq.simulate(ESEEM_depth=0.1, Tm=1e3)
plt.plot(x, data.real)
plt.xlabel("Interpulse Delay (ns)")
plt.ylabel("Signal Intensity (a.u.)")