import pyepr as epr
import matplotlib.pyplot as plt

pulse = epr.WURSTPulse(tp=128, flipangle=np.pi, init_freq=-0.2, final_freq=0.1,order=32)
pulse.plot(pad=100)