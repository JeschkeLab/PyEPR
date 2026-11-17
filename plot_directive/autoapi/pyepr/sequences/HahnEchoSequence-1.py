import pyepr as epr
import numpy as np

epr.set_waveform_precision(1)
seq = epr.HahnEchoSequence(
    B=12000, freq=34, reptime=1e3, averages=1, shots=10, tau=500)
seq.draw()