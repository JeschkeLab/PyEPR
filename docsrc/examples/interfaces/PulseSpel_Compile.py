"""
PulseSpel Compile 
===============================================

PyEPR can compile a Bruker PulseSpel pulse sequence from a PyEPR `sequence` object. 
This allows for easy integration of PyEPR with Bruker Xepr software, 
and allows for the use of PyEPR's powerful sequence building tools to create complex pulse sequences for EPR experiments.
"""

#%%
# Import packages
# --------------------------------------------
import pyepr as epr
from pyepr.hardware.Bruker_tools import write_pulsespel_file
import numpy as np
from pathlib import Path

#%%
# Build Sequence
# --------------------------------------------
# In this example we will build a simple 4 pulse Carr-Purcell sequence with a complete phase cycle. 

p90 = epr.RectPulse(tp=32,scale=0.1,flipangle=np.pi/2)
p180 = epr.RectPulse(tp=32,scale=0.2,flipangle=np.pi) 
seq = epr.CarrPurcellSequence(
    B=12000, freq=34, reptime=3e3, averages=1, shots=10,
    start=300, step=16, dim=512,n=4,
    pi2_pulse=p90, pi_pulse=p180)
print(seq)
seq.draw()


#%%
# Convert to PulseSpel and Save to File
# --------------------------------------------


def_file, exp_file, shapes = write_pulsespel_file(
    sequence=seq,
    d0 = 617, #ns
    AWG = True, # Use AWG or MPFU channels
)

out_dir = Path('_pulsespel_output')
out_dir.mkdir(exist_ok=True)
(out_dir / 'CP4.def').write_text(def_file)
(out_dir / 'CP4.exp').write_text(exp_file)
(out_dir / 'CP4.shp').write_text('\n'.join(shapes))
# sphinx_gallery_capture_repr_block = ()

#%%
# The generated files are:
# These are not as human-readable as commonly used PulseSpel files, but this isn't an issue. 
#
# .. tab-set::
#
#     .. tab-item:: Definitions (.def)
#
#         .. literalinclude:: /examples/interfaces/_pulsespel_output/CP4.def
#             :language: none
#             :class: scroll-40
#
#     .. tab-item:: Experiment (.exp)
#
#         .. literalinclude:: /examples/interfaces/_pulsespel_output/CP4.exp
#             :language: none
#             :class: scroll-40
#
#     .. tab-item:: Shapes (.shp)
#
#         .. literalinclude:: /examples/interfaces/_pulsespel_output/CP4.shp
#             :language: none
#             :class: scroll-40

#%%
# Limitations
# --------------------------------------------
#
# Currently known limitations of the PulseSpel compiler include:
#
# - Does not support RF pulses (only microwave pulses are supported)
# - Does not support transient detection (only echo detection is supported)
