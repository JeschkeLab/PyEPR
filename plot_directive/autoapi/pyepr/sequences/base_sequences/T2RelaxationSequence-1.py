import pyepr as epr
epr.set_waveform_precision(1)
seq = epr.T2RelaxationSequence(
    B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=300, step=40, dim=200)
seq.draw()