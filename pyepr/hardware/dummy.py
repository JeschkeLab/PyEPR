from pyepr.classes import  Interface, Parameter
from pyepr.dataset import  create_dataset_from_sequence
from pyepr.pulses import Pulse, RectPulse, ChirpPulse, HSPulse, Delay, Detection
from pyepr.sequences import *
from pyepr.fieldsweep_analysis import create_Nmodel
import yaml

import numpy as np
import deerlab as dl
import time
import logging
from dataclasses import dataclass
import os

rng = np.random.default_rng(12345)

hw_log = logging.getLogger('interface.Dummy')

def val_in_us(Param):
        if len(Param.axis) == 0:
            if Param.unit == "us":
                return Param.value
            elif Param.unit == "ns":
                return Param.value / 1e3
        elif len(Param.axis) == 1:
            if Param.unit == "us":
                return Param.tau1.value + Param.axis[0]['axis']
            elif Param.unit == "ns":
                return (Param.value + Param.axis[0]['axis']) / 1e3 

def val_in_ns(Param):
        if len(Param.axis) == 0:
            if Param.unit == "us":
                return Param.value * 1e3
            elif Param.unit == "ns":
                return Param.value 
        elif len(Param.axis) == 1:
            if Param.unit == "us":
                return (Param.tau1.value + Param.axis[0]['axis']) * 1e3
            elif Param.unit == "ns":
                return (Param.value + Param.axis[0]['axis']) 

def add_noise(data, noise_level):
    # Add noise to the data with a given noise level for data that could be either real or complex
    if np.isrealobj(data):
        noise = np.squeeze(rng.normal(0, noise_level, size=(*data.shape,1)).view(np.float64))

    else:
        noise = np.squeeze(rng.normal(0, noise_level, size=(*data.shape,2)).view(np.complex128))
    data = data + noise
    return data

def add_phaseshift(data, phase):
    data = data.astype(np.complex128) * np.exp(-1j*phase*np.pi)
    return data
    
@dataclass()
class dummySample:
    """
    This is a data class that holdes the information for a dummy sample allowing for the simulation of the experiment.

    Parameters:
    -----------
    name: str
        The name of the sample
    conc: float
        The spin concentration of the sample, this is used for both signal strength and DEER simulation
    Tm: float
        The Tm of the sample
    
    """

    name: str
    conc:float        

@dataclass()
class dummyResonator:
    fc:float
    Q:float
    nu1:float
    noise_level:float

    def mode(self,x):
        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y

        mode = lambda x: lorenz_fcn(x, self.fc, self.fc/self.Q)
        scale = self.nu1/mode(self.fc)

        return lorenz_fcn(x, self.fc, self.fc/self.Q) * scale

        

class dummyInterface(Interface):


    def __init__(self,config_file) -> None:
        """
        Parameters:
        -----------
        config_file: str or dict
            Either the path to the config file or a dictionary containing the config file. The config file should be a yaml file containing the following
            keys:
            - Spectrometer
                - Dummy: 
                    - speedup: float
                    - SNR: float
                    - ESEEM_depth: float
                    - noise_level: float
                    - phaseshift: float
                    - Sample:
                        - name: str
                        - conc: float
                        - Tm: float
                - Bridge: 
                    
        """
        
        if isinstance(config_file, dict):
            config = config_file
            self.config = config

        elif os.path.exists(config_file) == True:
            with open(config_file, mode='r') as file:
                config = yaml.safe_load(file)
                self.config = config
        else:
            raise FileNotFoundError("Config file not found")
        
        Dummy = config['Spectrometer']['Dummy']
        Bridge = config['Spectrometer']['Bridge']
        resonator_list = list(config['Resonators'].keys())
        self.state = False
        self.speedup = Dummy['speedup']
        self.pulses = {}
        self.start_time = 0
        self.SNR = Dummy['SNR']
        self.phaseshift = Dummy.get('phaseshift',0.1)
        self.ESEEM = Dummy.get('ESEEM_depth',0)


        # Create virtual Resonator
        key1 = resonator_list[0]
        fc = self.config['Resonators'][key1]['Center Freq']
        Q = self.config['Resonators'][key1]['Q']
        nu1 = self.config['Resonators'][key1].get('nu1',75)
        noise_level = Dummy.get('noise_level',0.02)
        
        self.dummyResonator = dummyResonator(fc,Q,nu1,noise_level)

        # Create virtual Sample
        if 'DummySample' in Dummy.keys():
            DummySample = Dummy['Sample']
        else:
            DummySample = {'name':'dummy','conc':20,'Tm':1.5}
        name = DummySample.get('name','dummy')
        conc = DummySample.get('conc',20)
        Tm = DummySample.get('Tm',1.5) #us

        self.dummySample = dummySample(name,conc)


        self.mode = self.dummyResonator.mode
        super().__init__(log=hw_log)

    def launch(self, sequence, savename: str, **kwargs):
        hw_log.info(f"Launching {sequence.name} sequence")
        self.state = True
        self.cur_exp:Sequence = sequence
        self.start_time = time.time()
        return super().launch(sequence, savename)
    
    def acquire_dataset(self,**kwargs):
        hw_log.debug("Acquiring dataset")

        if hasattr(self.cur_exp,'simulate'):
            axes, data = self.cur_exp.simulate()
        else:
            raise NotImplementedError("Simulation not implemented for this sequence")

        time_estimate = self.cur_exp._estimate_time()
        if self.speedup != np.inf:
            time_estimate /= self.speedup
            progress = (time.time() - self.start_time) / time_estimate
            if progress > 1:
                progress = 1
            elif progress < 0.01:
                progress = 0.01
            n_acq = int(progress * self.cur_exp.averages.value * self.cur_exp.shots.value * self.cur_exp.pcyc_dets.shape[0])
            if n_acq == 0:
                n_acq = 1
            SNR = ((self.dummySample.conc * 1e-3) / self.dummyResonator.noise_level) * np.sqrt(n_acq)

            data = add_noise(data, 1/SNR)
        else:
            progress = 1

        data = add_phaseshift(data,self.phaseshift)
        scan_num = self.cur_exp.averages.value
        dset = create_dataset_from_sequence(data,self.cur_exp)
        dset.attrs['nAvgs'] = int(scan_num*progress)
        
    
        return super().acquire_dataset(dset)
    
    def tune_rectpulse(self,*,tp, freq, B, reptime,**kwargs):

        rabi_freq = self.mode(freq)
        def Hz2length(x):
            return 1 / ((x/1000)*2)
        rabi_time = Hz2length(rabi_freq)
        if rabi_time > tp:
            p90 = tp
            p180 = tp*2
        else:
            p90 = rabi_time/tp
            p180 = p90*2

        self.pulses[f"p90_{tp}"] = RectPulse(tp=tp, freq=0, flipangle=np.pi/2, scale=p90)
        self.pulses[f"p180_{tp}"] = RectPulse(tp=tp, freq=0, flipangle=np.pi, scale=p180)

        return self.pulses[f"p90_{tp}"], self.pulses[f"p180_{tp}"]
    
    def tune_pulse(self, pulse, mode, freq, B , reptime, shots=400):
        hw_log.debug(f"Tuning {pulse.name} pulse")
        pulse.scale = Parameter('scale',0.5,unit=None,description='The amplitude of the pulse 0-1')
        hw_log.debug(f"Setting {pulse.name} pulse to {pulse.scale.value}")
        return pulse
            
    def isrunning(self) -> bool:
        current_time = time.time()
        runtime =  (self.cur_exp._estimate_time() / self.speedup)
        runtime = np.min([runtime, 5])
        if current_time - self.start_time > runtime:
            self.state = False
    
        return self.state
    
    def terminate(self) -> None:
        self.state = False
        hw_log.info("Terminating sequence")
        return super().terminate()
    

