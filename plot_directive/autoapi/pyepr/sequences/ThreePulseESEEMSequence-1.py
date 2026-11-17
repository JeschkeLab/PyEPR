import pyepr as epr
epr.set_waveform_precision(1)
seq = epr.ThreePulseESEEMSequence(
    B=12000, freq=34, reptime=1000, averages=1, shots=10,
    tau=300, start=500, step=40, dim=200)
seq.draw()