import pyepr as epr
epr.set_waveform_precision(1)
seq = epr.CarrPurcellSequence(B=12000, freq=34, reptime=1000, averages=1, shots=10, n=4)

seq.draw()