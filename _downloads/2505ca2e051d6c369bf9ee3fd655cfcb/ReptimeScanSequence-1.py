import pyepr as epr
epr.set_waveform_precision(1)
seq = epr.ReptimeScanSequence(
B=12000, freq=34, reptime=1000, averages=1,
shots=10, reptime_max=10000, start=20, dim=100)

seq.draw()