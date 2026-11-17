import pyepr as epr
import numpy as np
epr.set_waveform_precision(1)

seq = epr.TWTProfileSequence(
    B=12000, freq=34, reptime=3e3, averages=1, shots=10, dtp=2)
seq.draw()