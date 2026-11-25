"""
Uhrig Decoupling Sequence
===============================================

The Uhrig dynamical decoupling (UDD) sequence uses ``n`` refocusing pulses
placed non-uniformly over the total evolution time ``T``. The j-th pulse is
placed at

.. math::
    t_j = T \\sin^2\\left(\\frac{\\pi j}{2n + 2}\\right), \\quad j = 1, \\dots, n

and the echo forms at ``T``.
"""
# %%
import pyepr as epr
import numpy as np
import matplotlib.pyplot as plt

epr.set_waveform_precision(1)

# %%
# Define the total evolution time ``T`` as a virtual parameter, this is
# stepped but never passed to the spectrometer itself.

n = 5  # Number of refocusing pulses
T = epr.Parameter(
    "T", 1000, step=200, dim=50, unit="ns",
    description="Total evolution time", virtual=True)

# %%
# Build the sequence. Each refocusing pulse position is a fixed fraction of
# ``T``, so its step is generally not a multiple of the waveform precision.
# Rounding the step would accumulate errors along the axis, so instead each
# position is given as an explicit (non-uniform) axis, rounded point by point,
# and linked to the ``T`` axis.
#
# The same phase cycle as the Carr-Purcell sequence is used. Shifting the
# phase of a refocusing pulse by :math:`\phi` shifts the echo phase by
# :math:`2\phi`, which sets the detection signs.
#
# ======================================  ==============  ===============
# Pulse                                   Phases          Detection signs
# ======================================  ==============  ===============
# :math:`\pi/2`                           +x, -x          +, -
# :math:`\pi_1` to :math:`\pi_{n-2}`      +x, +y, -x, -y  +, -, +, -
# :math:`\pi_{n-1}`                       +x, -x          +, +
# :math:`\pi_n`                           +x              \+
# ======================================  ==============  ===============
#
# The phase cycles are combined multiplicatively, giving
# :math:`2 \times 4^{n-2} \times 2` steps in total.

seq = epr.Sequence(
    name="UhrigSequence", B=12127, freq=34, reptime=3e3, averages=1, shots=10)

seq.addPulse(epr.RectPulse(
    t=0, tp=16, freq=0, flipangle=np.pi/2,
    pcyc={"phases": [0, np.pi], "dets": [1, -1]}))

for j in range(1, n + 1):
    if j == n:
        pcyc = {"phases": [0], "dets": [1]}
    elif j == n - 1:
        pcyc = {"phases": [0, np.pi], "dets": [1, 1]}
    else:
        pcyc = {"phases": [0, np.pi/2, np.pi, -np.pi/2], "dets": [1, -1, 1, -1]}

    t_j = epr.Parameter(
        f"t{j}", np.round(T.get_axis() * np.sin(np.pi * j / (2 * n + 2))**2),
        unit="ns", link=T, description=f"Position of refocusing pulse {j}")
    seq.addPulse(epr.RectPulse(
        t=t_j, tp=32, freq=0, flipangle=np.pi, pcyc=pcyc))

seq.addPulse(epr.Detection(t=T, tp=512))
seq.evolution([T])

print(seq)
print(f"Number of phase cycle steps: {seq.pcyc_dets.shape[0]}")

# %%
# Plot the pulse positions as a fraction of ``T``, compared to a Carr-Purcell
# sequence with the same number of pulses.

udd = [seq.pulses[j].t.value / T.value for j in range(1, n + 1)]
cp = [(2 * j - 1) / (2 * n) for j in range(1, n + 1)]

plt.figure(figsize=(6, 2))
plt.eventplot([cp, udd], lineoffsets=[0, 1], colors=["C0", "C1"])
plt.yticks([0, 1], ["Carr-Purcell", "Uhrig"])
plt.xlim(0, 1)
plt.xlabel("t / T")
plt.tight_layout()

# %%
seq.draw()
