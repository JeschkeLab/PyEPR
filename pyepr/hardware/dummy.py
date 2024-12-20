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
    

class dummyInterface(Interface):


    def __init__(self,config_file) -> None:
        with open(config_file, mode='r') as file:
            config = yaml.safe_load(file)
            self.config = config
        
        Dummy = config['Spectrometer']['Dummy']
        Bridge = config['Spectrometer']['Bridge']
        resonator_list = list(config['Resonators'].keys())
        self.state = False
        self.speedup = Dummy['speedup']
        self.pulses = {}
        self.start_time = 0
        self.SNR = Dummy['SNR']
        if 'ESEEM_depth' in Dummy.keys():
            self.ESEEM = Dummy['ESEEM_depth']
        else:
            self.ESEEM = 0

        # Create virtual mode
        key1 = resonator_list[0]
        fc = self.config['Resonators'][key1]['Center Freq']
        Q = self.config['Resonators'][key1]['Q']
        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y

        mode = lambda x: lorenz_fcn(x, fc, fc/Q)
        x = np.linspace(Bridge['Min Freq'],Bridge['Max Freq'])
        scale = 75/mode(x).max()
        self.mode = lambda x: lorenz_fcn(x, fc, fc/Q) * scale
        super().__init__(log=hw_log)

    def launch(self, sequence, savename: str, **kwargs):
        hw_log.info(f"Launching {sequence.name} sequence")
        self.state = True
        self.cur_exp = sequence
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
            data = add_noise(data, 1/(self.SNR*progress))
        else:
            progress = 1
        scan_num = self.cur_exp.averages.value
        dset = create_dataset_from_sequence(data,self.cur_exp)
        dset.attrs['nAvgs'] = int(scan_num*progress)
        
    
        return super().acquire_dataset(dset)
    
    def tune_rectpulse(self,*,tp, LO, B, reptime,**kwargs):

        rabi_freq = self.mode(LO)
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
    
    def tune_pulse(self, pulse, mode, LO, B , reptime, shots=400):
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
    

