import pyepr as epr
import numpy as np
epr.set_waveform_precision(1)

seq = epr.Sequence(
    name="CustomHahnEcho", B=12000, freq=34, reptime=1e3,
    averages=1, shots=10)
seq.addPulse(epr.RectPulse(t=0, tp=16, freq=0, flipangle=np.pi/2))
seq.addPulse(epr.RectPulse(t=500, tp=32, freq=0, flipangle=np.pi))
seq.addPulse(epr.Detection(t=1000, tp=512))
seq.draw()