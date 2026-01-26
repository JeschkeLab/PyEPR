# Interfacing with Hardware

PyEPR has been designed to interface with a range of different EPR spectrometers, providing a hardware abstraction layer that allows users to write experiments in a hardware-agnostic manner. Currently, PyEPR supports the following hardware interfaces:
    - Bruker ElexSys AWG based spectrometers
    - ETH Zürich Matlab based homebuilt spectrometers

For homebuilt spectrometers, there can be additional packages required. Please refer to the respective hardware interface documentation for more details.

## Initializing Hardware Interfaces
The hardware interfaces can be initialized by importing the respective classes from the `pyepr.hardware` module and providing a configuration file that specifies the hardware settings.
```{eval-rst}
.. tab-set::

    .. tab-item:: Bruker AWG
        
        .. code-block:: python

            from pyepr.hardware.Bruker_AWG import BrukerAWG
            # Initialize the Bruker AWG interface with a configuration file
            interface = BrukerAWG('path/to/config_file.yaml')

    .. tab-item:: ETH Matlab
        
        .. code-block:: python

            from pyepr.hardware.ETH_awg import ETH_awg_interface
            # Initialize the ETH Matlab interface with a configuration file
            interface = ETH_awg_interface('path/to/config_file.yaml')


```

### Configuration Files

```{eval-rst}
.. tab-set::

    .. tab-item:: Bruker AWG
        
        .. literalinclude:: examples/config_files/BrukerElexSys_config.yaml
            :language: yaml
        

    .. tab-item:: ETH Matlab
        
        .. literalinclude:: examples/config_files/ETHmatlab_config.yaml
            :language: yaml

    .. tab-item:: Dummy

        .. literalinclude:: examples/config_files/Dummy_config.yaml
            :language: yaml
        
```

### First Tests
After initializing the hardware interface, it is recommended to run some basic tests to ensure that the connection to the spectrometer is working correctly. 



## Pulse Tuning

The hardware interface can be used to create and tune pulses for the spectrometer. 

### Creating rectangular pulse pairs
Often we just need a pair of rectangular pulses. These can be quickly created though an amplitude sweep. 

```python
p90, p180 = interface.tune_rectpulse(
                                     tp = 16,  # Pulse length of pi/2 pulse in ns
                                     B = 12200,  # Magnetic field in Gauss
                                     freq = 34.0,  # Microwave frequency in GHz
                                     reptime = 3000,  # Repetition time in us
                                     shots = 20,  # Number of shots per point
                                     )
```
This will create and return a pi/2 and pi rectangular pulse tuned to the specified field and frequency. The pulses can then be used directly in sequences, and have the correct scale (amplitude) set for the hardware.

### Tuning pre-defined pulse shapes
Often a specific pulse shape is required, and it is better to tune that pulse directly. This can be done using the `tune_pulse` method of the hardware interface. 
This is done as an amplitude sweep, either as a Hahn Echo (`amp_hahn`) or more commonly as a hole-buring recovery experiment (`amp_nut`). 
When a hole-burning recovery experiment is used, a pair of on-resonance rectangular pulses are first created using the `tune_rectpulse` method.

```python
from pyepr.pulses import GaussPulse
pulse = GaussPulse(tp=32, freq=0, flipangle=np.pi)  # Create a Gaussian pulse with length 32 ns
tuned_pulse = interface.tune_pulse(
                                     pulse,
                                     'amp_nut',
                                     B = 12200,  # Magnetic field in Gauss
                                     freq = 34.0,  # Microwave frequency in GHz
                                     reptime = 3000,  # Repetition time in us
                                     shots = 20,  # Number of shots per point
                                     )
```