import pyepr as epr
import matplotlib.pyplot as plt

pulse = epr.SincPulse(tp=128, flipangle=np.pi, freq=0, order=6)
pulse.plot(pad=100)