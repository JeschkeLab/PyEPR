import pyepr as epr
epr.set_waveform_precision(1)

seq = epr.FieldSweepSequence(
B=12000, Bwidth=100, freq=34, reptime=1e3, averages=1, shots=10)
seq.draw()