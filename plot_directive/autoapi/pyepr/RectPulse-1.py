import pyepr as epr
import matplotlib.pyplot as plt

pulse = epr.RectPulse(tp=32, flipangle=np.pi)
pulse.plot(pad=50)