# External imports
import numpy as np
import time
import re
import warnings
import requests
import pickle
import datetime
import logging

# PyEPR imports
from pyepr.classes import  Interface, Parameter
from pyepr.pulses import  Delay, Detection
from pyepr.sequences import Sequence, HahnEchoSequence


log = logging.getLogger("interface")
class PyEPRControlInterface(Interface):

    """
    Represents the interface for spectrometers built using PyEPR and connected with a local flask server.
    """


    def __init__(self,config_file_path:str):
        """

        Parameters
        -----------
        config_file_path : float
            The file path for the spectrometer configuration file.
        """
        
        
        super().__init__()
        self.IFgain_options = np.array([0, 20, 40])
        self.IFgain = 2

        self.config_file = config_file_path
        self.server = None

    @property
    def savefolder(self):
        return self._savefolder
    
    @savefolder.setter
    def savefolder(self, folder):
        self._savefolder = folder
        if hasattr(self,'server') and self.server is not None:
            response = requests.post(self.server + '/set_savefolder', json={'savefolder': folder})
            if 'error' in response.json():
                raise RuntimeError(response.json()['error'])
    
    def connect(self,ip='localhost',port=5000):
        ## Connect to the server

        self.ip = ip
        self.port = port
        self.server = f'http://{ip}:{port}'
        response = requests.post(self.server + '/connect', json={'config_file': self.config_file})

        if 'error' in response.json():
            raise RuntimeError(response.json()['error'])
        else:
            log.info(response.json()['message'])
        
        if self.savefolder is not None:
            self.savefolder = self.savefolder
        
    def disconnect(self):
        response = requests.post(self.server + '/disconnect')
        log.info(response.json()['message'])
        self.server = None
        return True

    def acquire_dataset(self,verbosity=0,sum_scans=True):
        
        for i in range(60):
            response = requests.get(self.server + '/get_data')
            if 'error' in response.json():
                if verbosity > 0:
                    log.warning(response.json()['error'])
                time.sleep(5)
                continue
            elif 'data' in response.json():
                break
            time.sleep(2)
        
        if 'data' in response.json():
            data = pickle.loads(response.json()['data'].encode('latin1')) #xr.datarray
            if sum_scans and 'Scan' in data.dims:
                return data.sum('Scan',keep_attrs = True)
            else:
                return data
        else:
            raise RuntimeError("No data returned from server")
    
    def terminate(self):
        response = requests.post(self.server + '/terminate')
        log.debug(response.json()['message'])
        return True
        

    def launch(self, sequence: Sequence , savename: str, IFgain=None, *args,**kwargs):

        # increase the detection length to a minimum 1024ns
        for pulse in sequence.pulses:
            if isinstance(pulse, Detection):
                if pulse.tp.value < 1024:
                    pulse.tp.value = 1024
        
        if (IFgain is None) or (IFgain is True):
            test_IF = True
            while test_IF:
                self.terminate()
                self._launch(sequence,savename,self.IFgain)
                scan_time = sequence._estimate_time() / sequence.averages.value
                check_1stScan = True
                while check_1stScan:
                    time.sleep(10)
                    dataset = self.acquire_dataset()
                    # if dataset.nAvgs == 0:
                    #     time.sleep(np.min([scan_time//10,2]))
                    #     continue

                    dig_level = dataset.attrs['diglevel']
                    pos_levels = dig_level * self.IFgain_options / self.IFgain_options[self.IFgain]
                    pos_levels[pos_levels > 0.85] = 0
                    if dig_level == 0:
                        continue
                    if (pos_levels[self.IFgain] > 0.85) or  (pos_levels[self.IFgain] < 0.03):
                        best_IFgain = np.argmax(pos_levels)
                    else:
                        best_IFgain = self.IFgain

                    if np.all(pos_levels==0):
                        log.critical('Saturation detected with IF gain 0. Please check the power levels.')
                        raise RuntimeError('Saturation detected with IF gain 0. Please check the power levels.')
                    elif (best_IFgain < self.IFgain) and (dataset.nAvgs == 0):
                        new_IFgain = np.max([self.IFgain - 1,best_IFgain])
                        log.info(f"IF gain changed from {self.IFgain} to {new_IFgain}")
                        self.IFgain = new_IFgain
                        check_1stScan = False
                    elif (best_IFgain != self.IFgain) and (dataset.nAvgs >= 1):
                        new_IFgain = np.min([self.IFgain +1,best_IFgain])
                        log.info(f"IF gain changed from {self.IFgain} to {new_IFgain}")
                        self.IFgain = new_IFgain
                        check_1stScan = False
                    elif dataset.nAvgs >=1:
                        log.debug(f"IF gain {self.IFgain} is optimal")
                        check_1stScan = False
                        test_IF = False
        elif (IFgain is not None) and isinstance(IFgain, [int, float]):
            
            self._launch(sequence,savename,IFgain)

        else:
            raise ValueError(f"IFgain must be of type [None, bool, int, float]. {IFgain} is not valid.")

    def _launch(self, sequence: Sequence , savename: str, IFgain=0,reset_cur_exp=True):

        timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
        self.savename = timestamp + savename + '.h5'
        if reset_cur_exp:
            self.cur_exp = sequence

        new_IFgain = self.IFgain_options[IFgain]
        if isinstance(new_IFgain, np.ndarray):
            log.warning(f"IF gain {new_IFgain} is not a valid option. Using {new_IFgain[-1]}")
            new_IFgain = new_IFgain[-1]
        
        log.debug(f"Launching sequence {savename} with IF gain {new_IFgain}")
        response = requests.post(self.server + '/launch', data= pickle.dumps({'seq':sequence,'savefile':savename,'IFgain':new_IFgain}))

        if 'error' in response.json():
            raise RuntimeError(response.json()['error'])
        pass
    def isrunning(self) -> bool:
        response = requests.get(self.server + '/isrunning')

        return response.json()['isrunning']

    def tune_rectpulse(self,*,tp, freq, B, reptime, shots=400):
        """Generates a rectangular pi and pi/2 pulse of the given length at 
        the given field position. This value is stored in the pulse cache. 

        Parameters
        ----------
        tp : float
            Pulse length of pi/2 pulse in ns
        freq : float
            Central frequency of this pulse in GHz
        B : float
            Magnetic B0 field position in Gauss
        reptime: float
            Shot repetion time in us.
        shots: int
            The number of shots

        Returns
        -------
        p90: RectPulse
            A tuned rectangular pi/2 pulse of length tp
        p180: RectPulse
            A tuned rectangular pi pulse of length tp
        """

        amp_tune =HahnEchoSequence(
            B=B, freq=freq, reptime=reptime, averages=1, shots=shots
        )

        scale = Parameter("scale",0,dim=45,step=0.02)
        amp_tune.pulses[0].tp.value = tp
        amp_tune.pulses[0].scale = scale
        amp_tune.pulses[1].tp.value = tp * 2
        amp_tune.pulses[1].scale = scale

        amp_tune.evolution([scale])
        
        self.launch(amp_tune, "autoDEER_amptune")

        while self.isrunning():
            time.sleep(10)
        dataset = self.acquire_dataset()
        dataset = dataset.epr.correctphase

        data = np.abs(dataset.data)
        scale = np.around(dataset.pulse0_scale[data.argmax()].data,2)
        log.debug(f"Optimal scale for {tp} ns pulse is {scale}")
        if scale > 0.9:
            raise RuntimeError("Not enough power avaliable.")
        
        if scale == 0:
            warnings.warn("Pulse tuned with a scale of zero!")
        p90 = amp_tune.pulses[0].copy(
            scale=scale, freq=amp_tune.freq)
        
        p180 = amp_tune.pulses[1].copy(
            scale=scale, freq=amp_tune.freq)

        return p90, p180

    
    def tune_pulse(self, pulse, mode, freq, B , reptime, shots=400):
        """Tunes a single pulse a range of methods.

        Parameters
        ----------
        pulse : Pulse
            The Pulse object in need of tuning.
        mode : str
            The method to be used.
        freq : float
            The local oscilator frequency in GHz
        B : float
            Magnetic B0 field position in Gauss
        reptime : us
            Shot repetion time in us.
        shots: int
            The number of shots

        Returns
        -------
        Tunned Pulse: Pulse
            The returned pulse object that is now tunned.
        """
        # Check pulse is a pulse
        if type(pulse) == Delay:
            pass
        if type(pulse) == Detection:
            pass
        
        # Get absolute central frequency
        if hasattr(pulse,"freq"):
            c_frq = pulse.freq.value + freq
        elif hasattr(pulse, "init_freq") & hasattr(pulse, "BW"):
            c_frq = pulse.init_freq.value + 0.5*pulse.BW.value + freq
        elif hasattr(pulse, "final_freq") & hasattr(pulse, "BW"):
            c_frq = pulse.final_freq.value - 0.5*pulse.BW.value + freq
        elif hasattr(pulse, "init_freq") & hasattr(pulse, "final_freq"):
            c_frq = 0.5*(pulse.final_freq.value + pulse.final_freq.value) + freq

        # Find rect pulses
        if mode == "amp_hahn":
            if pulse.flipangle.value == np.pi:
                tp = pulse.tp.value / 2
            elif pulse.flipangle.value == np.pi/2:
                tp = pulse.tp.value

            pi2_pulse, pi_pulse = self.tune_rectpulse(tp=tp, B=B, freq=c_frq, reptime=reptime)
            amp_tune =HahnEchoSequence(
                B=B, freq=freq, 
                reptime=reptime, averages=1, shots=shots,
                pi2_pulse = pulse, pi_pulse=pi_pulse
            )

            scale = Parameter('scale',0,unit=None,step=0.02, dim=45, description='The amplitude of the pulse 0-1')
            amp_tune.pulses[0].scale = scale

            amp_tune.evolution([scale])

            self.launch(amp_tune, "autoDEER_amptune")

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            dataset = dataset.epr.correctphase
            data = np.abs(dataset.data)

            new_amp = np.around(dataset.pulse0_scale[dataset.data.argmax()].data,2)
            if new_amp > 0.9:
                raise RuntimeError("Not enough power avaliable.")
        
            if new_amp == 0:
                warnings.warn("Pulse tuned with a scale of zero!")

            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
            return pulse

        elif mode == "amp_nut":
            pi2_pulse, pi_pulse = self.tune_rectpulse(tp=12, B=B, freq=c_frq, reptime=reptime)
            nut_tune = Sequence(
                name="nut_tune", B=(B/freq*c_frq), freq=freq, reptime=reptime,
                averages=1,shots=shots
            )
            nut_tune.addPulse(pulse.copy(
                t=0, pcyc={"phases":[0],"dets":[1]}, scale=0))
            nut_tune.addPulse(
                pi2_pulse.copy(t=2e3,
                               pcyc={"phases":[0, np.pi],"dets":[1, -1]},
                               freq=c_frq-freq))
            nut_tune.addPulse(
                pi_pulse.copy(t=2.5e3, pcyc={"phases":[0],"dets":[1]},
                              freq=c_frq-freq))
            nut_tune.addPulse(Detection(t=3e3, tp=512, freq=c_frq-freq))

            scale = Parameter('scale',0,unit=None,step=0.02, dim=45, description='The amplitude of the pulse 0-1')
            nut_tune.pulses[0].scale = scale
            nut_tune.evolution([scale])


            # nut_tune.addPulsesProg(
            #     pulses=[0],
            #     variables=["scale"],
            #     axis_id = 0,
            #     axis= np.arange(0,0.9,0.02)
            # )
            self.launch(nut_tune, "autoDEER_amptune")

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            dataset = dataset.epr.correctphase
            data = dataset.data
            axis = dataset.pulse0_scale
            # data = correctphase(dataset.data)
            if data[0] < 0:
                data *= -1

            if np.isclose(pulse.flipangle.value, np.pi):
                new_amp = np.around(axis[data.argmin()].data,2)
            elif np.isclose(pulse.flipangle.value, np.pi/2):
                sign_changes = np.diff(np.sign(np.real(data)))
                new_amp = np.around(axis[np.nonzero(sign_changes)[0][0]].data,2)
            else:
                raise RuntimeError("Target pulse can only have a flip angle of either: ",
                                "pi or pi/2.")
            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
        
            return pulse
    
    def tune(self,*, sequence=None, mode="amp_hahn", freq=None, gyro=None):

        if mode == "rect_tune":
            if freq is None:
                raise ValueError("freq must be given for rect_tune")
            if gyro is None:
                raise ValueError("gyro must be give")
            elif gyro >1:
                raise ValueError("Gyromagnetic ratio must be give in GHz/G")
            
            amp_tune =HahnEchoSequence(
                B=freq/gyro, freq=freq, reptime=2e3, averages=1, shots=400
            )
            tp = 12
            amp_tune.pulses[0].tp.value = tp
            amp_tune.pulses[0].scale.value = 0
            amp_tune.pulses[1].tp.value = tp*2
            amp_tune.pulses[1].scale.value = 0
            
            amp_tune.addPulsesProg(
                pulses=[0,1],
                variables=['scale','scale'],
                axis_id=0,
                axis=np.arange(0,0.9,0.02),
            )

            self.launch(amp_tune, "autoDEER_amptune", IFgain=0)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            scale = np.around(dataset.pulse0_scale[dataset.data.argmax()].data,2)
            if scale > 0.9:
                raise RuntimeError("Not enough power avaliable.")
            
            self.pulses[f"p90_{tp}"] = amp_tune.pulses[0].copy(
                scale=scale, freq=amp_tune.freq)
            self.pulses[f"p180_{tp*2}"] = amp_tune.pulses[1].copy(
                scale=scale, freq=amp_tune.freq)
        
        elif mode == "amp_hahn":
            for pulse in sequence.pulses:
                if type(pulse) == Delay:
                    continue
                if type(pulse) == Detection:
                    continue

                all_pulses = list(self.pulses.keys())
                pulse_matches = []
                for pulse_name in all_pulses:
                    if not re.match(r"^p180_",pulse_name):
                        continue
                    if not np.abs((self.pulses[pulse_name].freq.value + self.pulses[pulse_name].freq.value) - (sequence.freq.value + pulse.freq.value)) < 0.01:
                        continue
                    pulse_matches.append(pulse_name)
                    
                    ps_length_best =1e6
                for pulse_name in pulse_matches:
                    ps_length = int(re.search(r"p180_(\d+)",pulse_name).groups()[0])
                    if ps_length < ps_length_best:
                        ps_length_best = ps_length
                
                pi_pulse = self.pulses[f"p180_{ps_length_best}"]
                

                amp_tune =HahnEchoSequence(
                    B=sequence.B.value, freq=sequence.freq.value, 
                    reptime=sequence.reptime.value, averages=1, shots=400,
                    pi2_pulse = pulse, pi_pulse=pi_pulse
                )
        
                amp_tune.pulses[0].scale.value = 0

                amp_tune.addPulsesProg(
                    pulses=[0],
                    variables=['scale'],
                    axis_id=0,
                    axis=np.arange(0,0.9,0.02),
                )

                self.launch(amp_tune, "autoDEER_amptune")

                while self.isrunning():
                    time.sleep(10)
                dataset = self.acquire_dataset()
                scale = np.around(dataset.pulse0_scale[dataset.data.argmax()].data,2)
                pulse.scale.value = scale

            return sequence
  
