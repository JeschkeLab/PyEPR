"""
Linking a pulse's phase cycle to another pulse
===============================================

When several pulses in a sequence are phase cycled, PyEPR combines them
multiplicatively by default: each pulse's phase cycle is stepped through
independently, giving the full outer product of every combination.

Sometimes two pulses should instead be stepped *in phase* with one
another, rather than being combined with all other phase-cycled pulses.
This example shows how to link one pulse's phase cycle to another's so
that they are stepped through together.
"""
# %%
import numpy as np
import pyepr as epr

# %%
# Create a pulse with its own two-step phase cycle.
exc_pulse = epr.RectPulse(
    tp=16, freq=0, flipangle=np.pi / 2, t=0,
    pcyc={"phases": [0, np.pi], "dets": [1, -1]})

# %%
# Link a second pulse's phase cycle to the first by passing the pulse
# itself as the ``pcyc`` argument. This copies ``exc_pulse``'s phases and
# detection signs, and links the two pulses together.
ref_pulse = epr.RectPulse(
    tp=32, freq=0, flipangle=np.pi, t=100, pcyc=exc_pulse)

det_event = epr.Detection(tp=32, t=200)

# %%
# Build a sequence and add the pulses.
seq = epr.Sequence(
    name="example", B=12127, freq=34, reptime=3e3, averages=1, shots=10)
seq.addPulse([exc_pulse, ref_pulse, det_event])

print(seq)

# %%
# Because ``ref_pulse`` is linked to ``exc_pulse``, the two pulses step
# through their phases together, giving only 2 phase-cycle shots instead
# of the 2 x 2 = 4 that would result from two independent 2-step phase
# cycles.
print(f"Number of phase cycle shots: {seq.pcyc_dets.shape[0]}")
