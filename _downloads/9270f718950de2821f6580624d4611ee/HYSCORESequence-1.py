import pyepr as epr
epr.set_waveform_precision(1)

seq = epr.HYSCORESequence(
    B=12000, freq=34, reptime=1000, averages=1, shots=10,
    tau=300, start=300, step=16, dim=128)
seq.draw()