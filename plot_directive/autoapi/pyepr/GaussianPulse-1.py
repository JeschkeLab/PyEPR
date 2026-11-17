import pyepr as epr
import matplotlib.pyplot as plt

GuassPulse = epr.GaussianPulse(tp=32, flipangle=np.pi)
fig = GuassPulse.plot(pad=50)
axs = fig.get_axes()
axs[0].annotate('', xy=(-GuassPulse.FWHM.value/2, 0.5), xytext=(GuassPulse.FWHM.value/2, 0.5), color='C1',arrowprops=dict(arrowstyle='<|-|>', color='C1',lw=2))
axs[0].annotate('FWHM', xy=(0, 0.51), xytext=(0, 0.51), color='C1', ha='center', fontsize=12)