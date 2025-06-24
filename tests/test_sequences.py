from pyepr.sequences import FieldSweepSequence, CarrPurcellSequence, T2RelaxationSequence
import pytest
import numpy as np

def test_fieldsweep_simulate():
    
    N_gyro = 0.002803632236095
    freq = 34
    seq = FieldSweepSequence(
        B=freq/N_gyro,
        freq=freq,
        Bwidth=200,
        reptime=3e3,
        averages=1,
        shots=10

    )

    axes, data = seq.simulate()


def test_CarrPurcellSequence():
    N_gyro = 0.002803632236095
    freq = 34

    seq = CarrPurcellSequence(
        B=freq/N_gyro,
        freq=freq,
        reptime=3e3,
        averages=1,
        shots=10,
        n=2,
        dim=10,
        step=500
    )

    axes, data = seq.simulate()
    assert np.allclose(axes,np.arange(300,5000,step=500))
    assert np.allclose(data,[0.96619141-0.15302968j, 0.86849221-0.13755565j,0.72573584-0.11494526j, 0.5677971 -0.08993023j,0.4176877 -0.06615523j, 0.28976763-0.04589468j,0.19001049-0.03009471j, 0.11798591-0.01868713j,0.06948146-0.01100478j, 0.03885595-0.00615418j])

def test_T2RelaxationSequence():
    N_gyro = 0.002803632236095
    freq = 34

    seq = T2RelaxationSequence(
        B=freq/N_gyro,
        freq=freq,
        reptime=3e3,
        averages=1,
        shots=10,
        dim=10,
        step=500
    )

    axes, data = seq.simulate()
    assert np.allclose(axes,np.arange(500,5100,step=500))
    assert np.allclose(data.real,[6.39143954e-01, 2.90680188e-01, 1.31219296e-01, 3.81219609e-02,1.16796898e-02, 2.39331041e-03, 5.31640969e-04, 8.06720436e-05,1.34965908e-05, 1.56414588e-06])