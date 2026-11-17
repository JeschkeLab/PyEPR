import pyepr as epr
import matplotlib.pyplot as plt

pulse = epr.HSPulse(tp=128, flipangle=np.pi, init_freq=-0.2, final_freq=0.1, order1=1, order2=6, beta=20)
pulse.plot(pad=100)