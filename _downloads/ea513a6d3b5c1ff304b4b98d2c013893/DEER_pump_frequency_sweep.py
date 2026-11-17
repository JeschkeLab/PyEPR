"""
Testing different pump-offsets in Four-Pulse DEER  
==================================================
"""

import pyepr as epr
import numpy as np

obs_freq = 0 # Relative to the center frequency
pump_freqs = [-0.12,-0.10,-0.08,-0.06]

exc_pulse = epr.RectPulse(tp=32,scale=0.1,flipangle=np.pi/2,freq=obs_freq)
ref_pulse = epr.RectPulse(tp=32,scale=0.2,flipangle=np.pi,freq=obs_freq)
pump_pulse = epr.RectPulse(tp=32,scale=0.2,flipangle=np.pi,freq=pump_freqs[0])

seq = epr.FourPulseDEERSequence(
    B=12000, freq=34, reptime=3e3, averages=50, shots=25,
    start=300, step=16, dim=512,
    tau1=0.4,tau2=4,
    pi2_pulse=exc_pulse, pi_pulse=ref_pulse, pump_pulse=pump_pulse)

print(seq)
seq.draw()


#%%
# Connect to the spectrometer
#----------------------------


#%%

for pump_freq in pump_freqs:
    print(f"Running DEER with pump offset frequency {pump_freq*1e3} MHz")
    seq.pulses[2].freq.value = pump_freq

    # interface.launch()