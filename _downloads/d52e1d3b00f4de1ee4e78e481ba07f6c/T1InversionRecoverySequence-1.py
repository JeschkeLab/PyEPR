import pyepr as epr
epr.set_waveform_precision(1)
seq = epr.T1InversionRecoverySequence(
    B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=10000, step=2500, dim=500)
seq.draw()