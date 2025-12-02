# Pulse Sequencer

PyEPR provides an intuitive object-oriented pulse programmer allowing the user to design pulsesequences in a hardware-agnostic manner. Additionally, several common EPR experiments are pre-defined and can be easily instantiated and modified.

PyEPR uses ns, GHz and G as the default time, frequency and field units. Very occasionally, other units such as µs or MHz are used, in which case it will be explicitly mentioned.

## Seqeunce Construction

First we must import the pyepr package, and set our waveform precision. This waveform precision determines the time resolution of the pulse length and positons. If values are given that are not multiples of the waveform precision, they will be rounded.

```python
import pyepr as epr
epr.set_waveform_precision(2)  # Set waveform precision to 2 ns
```

Next, we want to create a sequence object. Here we are aiming to create a simple Hahn Echo sequence.
It is normally recommened, to set the field `B` and frequency `freq` as external variables. 
```python
B = 12200  # Magnetic field in Gauss
freq = 34.0  # Microwave frequency in GHz
seq = epr.Sequence(name='Hahn Echo Sequence',
                   B = B,
                   freq = freq,
                   reptime = 3e3,# Repetition time in us
                   averages = 1, 
                   shots = 20, # Number of shots per point
                   ) 
```
Now we need to define some pulses that can be used in our sequence. Here we create a 90 degree and a 180 degree rectangular pulse.
These pulses will eventually need a scale (amplitude), before the sequence can be run on hardware.
A Detection window is also created
```python
p90 = epr.RectPulse(tp=16,
                    freq=0, # Frequency offset in MHz, w.r.t the sequence frequency,
                    flipangle=np.pi/2, # Flip angle in degrees
                    pcyc = {"phases":[0, np.pi], "dets":[1,-1]}
                    )
p180 = epr.RectPulse(tp=32,
                    freq=0, # Frequency offset in MHz, w.r.t the sequence frequency,
                    flipangle=np.pi, # Flip angle in degrees
                    )
det = epr.Detetction(tp=32,
                    freq=0, # Frequency offset in MHz, w.r.t the sequence frequency,
                    )
```                  
We now need a time axis for our sequence and to add them to the sequence object.
When a pulse is copied into the sequence using the `add_pulse` method, parameters can be modified allowing the same pulse can be used multiple times with different timings or amplitudes.
```python
t = epr.Parameter(name='Interpulse Delay',
                  value=400, # Initial interpulse delay in ns
                  step=8,  # Step size in ns
                  dim=1024  # Number of points,
                  unit='ns'  # Unit of the parameter
                  description='Interpulse delay between the pi/2 and pi pulse'
                  )

# Adding the pulses to the sequence
seq.add_pulse(p90.copy(t=0))
seq.add_pulse(p180.copy(t=t))
seq.add_pulse(det.copy(t=2*t))  

# Defining the evolution
seq.evolution([t])
```
### Advanced Sequences
The sequence class is capable of more advanced features such as:
- Linked axes
- Reduced axes
- Non-linear axes

## Default Experiments
For convienience, several common EPR experiments are pre-defined and can be easily instantiated and modified.
A list of currently implemented experiments can be found in the [API documentation](API_docs.rst).
Here we show how to create a simple Hahn Echo Relaxation experiment.
```python
HE_Seq = epr.HahnEchoRelaxationSequence(
    B = B,
    freq = freq,
    reptime = 3e3, # Repetition time in us
    averages = 1,
    shots = 20, # Number of shots per point
    start = 400, # Initial interpulse delay in ns
    step = 8,  # Step size in ns
    dim = 1024  # Number of points
    pi2_pulse = p90,  # The 90 degree pulse
    pi_pulse = p180  # The 180 degree pulse
)
```

## Advanced Pulses
More complex pulse shapes can be created using the pulse classes provided in the `pyepr.pulses` module. A list of currently implemented pulse shapes, and their necessary inputs can be found in the [API documentation](API_docs.rst).

Here we show how to create a simple chirp pulse
```python
from pyepr.pulses import ChirpPulse

pulse = ChirpPulse(
    tp = 128,  # Pulse length in ns
    init_freq = -0.1,  # Frequency offset in GHz, w.r.t the sequence frequency,
    final_freq = 0.1,  # Frequency offset in GHz, w.r.t the sequence frequency,
    flipangle = np.pi,  # Flip angle in radians
)
```
