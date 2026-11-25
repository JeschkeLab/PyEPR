"""
Advanced Sequence Axes
===============================================


"""
# %%
import pyepr as epr
import numpy as np

# %%
# Parameter Arithmetic
# --------------------
#
# Parameters can be added to, subtracted from and multiplied by numbers or
# other parameters. Any dynamic axes are carried through to the new parameter.

tau1 = epr.Parameter("tau1", 400, unit="ns", description="Interpulse delay 1")
tau2 = epr.Parameter("tau2", 2000, unit="ns", description="Interpulse delay 2")
t = epr.Parameter(
    "t", 200, step=8, dim=250, unit="ns", description="Pump pulse position")

pump_t = tau1 + t
ref2_t = 2 * tau1 + tau2

print(f"Pump pulse position: {pump_t.get_axis()[:5]} ... ns")
print(f"Second refocusing pulse position: {ref2_t.value} ns")

# %%
# Linked Axes
# -----------
#
# In experiments such as the resonator profile sequence two parameters need to
# be stepped together, we can do this by linking the two parameters together.
# Parameter Arithmetic automatically creates a new parameter that is linked to
# the original parameter.

fc = 34.0 # GHz
B = 12127.0 # Gauss
gyro = fc/B # GHz/Gauss

freq = epr.Parameter(
    "freq", fc, start=-0.3, step=0.02, dim=30, unit="GHz",
    description="Frequency")
B_field = epr.Parameter(
    "B", B, start=-0.3/gyro, step=0.02/gyro, dim=30, unit="Gauss",
    link=freq, description="B0 Field")

seq = epr.Sequence(
    name="linked_axes", B=B_field, freq=fc, reptime=3e3, averages=1, shots=10)
seq.freq = freq
seq.addPulse([
    epr.RectPulse(t=0, tp=16, freq=0, flipangle=np.pi/2),
    epr.RectPulse(t=500, tp=32, freq=0, flipangle=np.pi),
    epr.Detection(t=1000, tp=512)])
seq.evolution([freq])

# Both ``freq`` and ``B`` are now stepped along the same axis
print(f"Variables: {seq.progTable['Variable']}, axis IDs: {seq.progTable['axID']}")
print(seq)

# %%
# Reduced Axes
# ------------
# We can also tag one moving parameter to be automatically reduced to a single
# value, for example for tau1 averaging in 4p-DEER.

tau1 = epr.Parameter(
    "tau1", 400, step=16, dim=8, unit="ns", description="Interpulse delay 1")
tau2 = epr.Parameter("tau2", 2000, unit="ns", description="Interpulse delay 2")
t = epr.Parameter(
    "t", 200, step=8, dim=250, unit="ns", description="Pump pulse position")

seq = epr.Sequence(
    name="4pDEER_tau_avg", B=B, freq=fc, reptime=3e3, averages=1, shots=10)
seq.addPulse([
    epr.RectPulse(t=0, tp=16, freq=0, flipangle=np.pi/2),
    epr.RectPulse(t=tau1, tp=32, freq=0, flipangle=np.pi),
    epr.RectPulse(t=tau1 + t, tp=32, freq=-0.1, flipangle=np.pi),
    epr.RectPulse(t=2*tau1 + tau2, tp=32, freq=0, flipangle=np.pi),
    epr.Detection(t=2*tau1 + 2*tau2, tp=256)])
seq.evolution([t, tau1], reduce=[tau1])

print(f"Reduced: {seq.progTable['reduce']}")
print(seq)

# %%
# Non-linear Axes
# ----------------
# We can also define a parameter that is non-linear for non-uniform sampling
# for example.
#
# Non linear axes are not supported on all spectrometer interfaces, for example
# Bruker XeprAPI based systems.

tau = epr.Parameter(
    "tau", np.round(np.geomspace(200, 20000, 40)), unit="ns",
    description="Interpulse delay")

seq = epr.Sequence(
    name="T2_NUS", B=B, freq=fc, reptime=3e3, averages=1, shots=10)
seq.addPulse([
    epr.RectPulse(t=0, tp=16, freq=0, flipangle=np.pi/2),
    epr.RectPulse(t=tau, tp=32, freq=0, flipangle=np.pi),
    epr.Detection(t=2*tau, tp=256)])
seq.evolution([tau])

print(f"tau axis: {tau.get_axis()[:5]} ... ns")
print(seq)