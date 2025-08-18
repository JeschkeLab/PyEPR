from pyepr.classes import Interface
from pyepr.hardware.XeprAPI_link import XeprAPILink
from pyepr.hardware.Bruker_tools import PulseSpel, run_general, build_unique_progtable,write_pulsespel_file,step_parameters, get_specjet_data
from pyepr.dataset import *
from pyepr.utils import transpose_list_of_dicts, transpose_dict_of_list,round_step
from pyepr.pulses import RectPulse, Delay
from pyepr.sequences import HahnEchoSequence
import tempfile
import time
from scipy.optimize import minimize_scalar, curve_fit
import numpy as np
import threading
import concurrent.futures
from PyQt6.QtCore import QThreadPool
from autodeer.gui import Worker
import os
from pathlib import Path
import datetime
from deerlab import correctphase
import matplotlib.pyplot as plt
import logging
import warnings
# ===================


# =============================================================================
hw_log = logging.getLogger('interface.Xepr')


class BrukerAWG(Interface):
    """
    Represents the interface for connecting to AWG based Bruker ELEXSYS-II 
    Spectrometers.
    """

    def __init__(self, config_file:str) -> None:
        """An interface for connecting to AWG based Bruker ELEXSYS-II 
        Spectrometers.

        Getting Started
        ------------------
        Before a connection can be made an appropriate configuration file first
        needs to be written. 

        1. Open Xepr
        2. Processing -> XeprAPI -> Enable XeprAPI
        3. `BrukerAWG.connect()`

        Parameters
        ----------
        config_file : str
            The path to a YAML configuration file.
        
        Attributes
        ----------
        bg_thread: None or threading.Thread
            If a background thread is needed, it is stored here.

        """

        self.api = XeprAPILink(config_file)
        self.spec_config = self.api.config["Spectrometer"]
        self.bridge_config = self.api.spec_config["Bridge"]
        
        self.temp_dir = tempfile.mkdtemp("autoDEER")

        self.d0 = self.bridge_config["d0"]
        
        self.bg_thread = None
        self.bg_data = None
        self.cur_exp = None
        self.tuning = False
        self.pool = QThreadPool()
        self.savename = ''
        self.savefolder = str(Path.home())
        self.setup_flag=False

        super().__init__()


    def connect(self, d0=None) -> None:

        self.api.connect()

        time.sleep(1)
        self.setup(d0)
        
        return super().connect()

    def setup(self,d0=None):

        self.api.hidden['BrPlsMode'].value = True
        self.api.hidden['OpMode'].value = 'Operate'
        # self.api.hidden['RefArm'].value = 'On'

        # TODO add detection of VAMP-III model
        video_bw = self.bridge_config.get('Video BW')
        self.api.cur_exp['VideoBW'].value = 20

        self.time_base = 1/(video_bw*2e-3)
        self.api.hidden['specJet.TimeBase'].value = self.time_base


        if d0 is None:
            self.calc_d0()
        else:
            self.d0 = d0
        self.api.hidden['Detection'].value = 'Signal'
        self.setup_flag = True
        pass

    def acquire_dataset(self, **kwargs):
        if self.bg_data is None:
            if not self.isrunning():
                if (self.savename is not None) and (self.savename != ''):
                    full_path = os.path.join(self.savefolder,self.savename)
                    if Path(full_path).is_file():
                        self.savename = self.savename+"_1"
                    self.api.xepr_save(os.path.join(self.savefolder,self.savename))
                data = self.api.acquire_dataset(self.cur_exp)
            else:
                data = self.api.acquire_scan(self.cur_exp,**kwargs)
        else:

            data = create_dataset_from_sequence(self.bg_data, self.cur_exp)
        return super().acquire_dataset(data)


    def _launch_complex_thread(self,sequence,axID=1,tune=True):
    
        uProgTable = build_unique_progtable(sequence)
        reduced_seq = sequence.copy()
        reduced_seq.progTable = transpose_list_of_dicts([transpose_dict_of_list(sequence.progTable)[0]])

        uProgTable_py = uProgTable[axID]
        axis = uProgTable_py['axis']
        reduced_seq.averages.value = 1
        py_ax_dim = uProgTable_py['axis']['dim']

        
        self.bg_data = np.zeros((uProgTable[0]['axis']['dim'],py_ax_dim),dtype=np.complex128)
        
        print("Initial PulseSpel Launch")
        self.launch(reduced_seq,savename='test',tune=tune, update_pulsespel=True, start=False,reset_bg_data=False,reset_cur_exp=False)

        self.terminate()

        variables = uProgTable_py['variables']
        # print("Creating Thread")
        # # thread = threading.Thread(target=step_parameters,args=[self,reduced_seq,py_ax_dim,variables])
        # # self.bg_thread = self.pool.submit(step_parameters, self,reduced_seq,py_ax_dim,variables)
        # self.bg_thread = Worker(step_parameters, self,reduced_seq,py_ax_dim,variables)
        # self.pool.start(self.bg_thread)
        # print("Started Thread")

        step_parameters(self,reduced_seq,py_ax_dim,variables)
        # thread.start()

        pass
        
    def launch(self, sequence: Sequence, savename: str, start=True, tune=True,
            MPFU_overwrite=None,update_pulsespel=True, reset_bg_data = True,
            reset_cur_exp=True, IFgain=None, **kwargs):
        
        sequence.shift_detfreq_to_zero()

        if self.isrunning():
            self.terminate(now=True)
            time.sleep(4)
                    
        if reset_bg_data:
            self.bg_data = None

        if reset_cur_exp:
            timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
            if Path(os.path.join(self.savefolder,timestamp+savename+'.DSC')).is_file():
                timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M%S_')
            self.savename = timestamp+savename
            self.cur_exp = sequence
                    
        # First check if the sequence is pulsespel compatible
    
        if not test_if_MPFU_compatability(sequence):
            print("Launching complex sequence")
            self._launch_complex_thread(sequence,1,tune)
            return None
                
        
        pg = sequence.pulses[-1].tp.value
        pg = round_step(pg/2,self.bridge_config['Pulse dt'])
        d0 = self.d0-pg

        if d0 <0:
            d0=0
        d0 = round_step(d0,self.bridge_config['Pulse dt'])
        
        if update_pulsespel:

            MaxGate = np.min([40,sequence.reptime.value*self.bridge_config['DutyCycle']*1e-2])
            def_text, exp_text, pulse_shapes = write_pulsespel_file(sequence,d0,AWG=True,MaxGate=MaxGate)
            
            verbMsgParam = self.api.cur_exp.getParam('*ftEPR.PlsSPELVerbMsg')
            plsSPELCmdParam = self.api.cur_exp.getParam('*ftEPR.PlsSPELCmd')
            self.api.XeprCmds.aqPgSelectBuf(1)
            self.api.XeprCmds.aqPgSelectBuf(2)
            self.api.XeprCmds.aqPgSelectBuf(3)

            # Merge the pulse shapes into one string
            pulse_shapes = '\n'.join(pulse_shapes)

            self.api.cur_exp.getParam('*ftEpr.PlsSPELShpTxt').value = pulse_shapes

            plsSPELCmdParam.value=9 # Processing shapes
            time.sleep(2)
            self.api.XeprCmds.aqPgSelectBuf(1)

            self.api.cur_exp.getParam('*ftEpr.PlsSPELGlbTxt').value = def_text
            self.api.XeprCmds.aqPgShowDef()

            plsSPELCmdParam.value=3 # Compile variables
            time.sleep(2)
            # while not "The variable values are set up" in verbMsgParam.value:
            #     time.sleep(0.1)

            self.api.XeprCmds.aqPgSelectBuf(2)

            self.api.cur_exp.getParam('*ftEpr.PlsSPELPrgTxt').value = exp_text
            plsSPELCmdParam.value=5 # Compile experiments
            time.sleep(4)
            # while not "Second pass ended" in verbMsgParam.value:
            #     time.sleep(0.1)


                
        self.api.set_field(sequence.B.value)
        self.api.set_freq(sequence.freq.value)
        self.api.hidden['specJet.TimeBase'].value = self.time_base

        if IFgain is not None:
            self.api.set_video_gain(IFgain)
        self.IFgain = self.api.get_video_gain()
        
        self.api.set_PhaseCycle(self.bridge_config.get('On Board Pcyc',False))
        if 'B' in sequence.progTable['Variable']:
            idx = sequence.progTable['Variable'].index('B')
            B_axis = sequence.progTable['axis'][idx]
            self.api.set_sweep_width(B_axis.max()-B_axis.min())
            self.api.set_PhaseCycle(False)
            # self.api.set_field(sequence.B.value + B_axis.min())
            self.api.set_field(sequence.B.value)
            self.api.cur_exp['fieldCtrl.AtLeftBorder'].value = True
            time.sleep(5)
        

        self.api.set_ReplaceMode(False)
        self.api.set_Acquisition_mode(1)
        pg = sequence.pulses[-1].tp.value
        pg = round_step(pg/2,self.bridge_config['Pulse dt'])
        d0 = self.d0-pg

        if d0 <0:
            d0=0
        d0 = round_step(d0,self.bridge_config['Pulse dt'])
        time.sleep(1.0)
        self.api.set_PulseSpel_var('d0', d0)
        self.api.set_PulseSpel_var('pg', pg*2)
        self.api.set_PulseSpel_experiment('auto')
        self.api.set_PulseSpel_phase_cycling('auto')

        if start:
        # if start and not ('B' in sequence.progTable['Variable']):
        #     self.api.run_exp()
        # elif start:
        #     self.api.set_PulseSpel_var('m', 1)
        #     self.api.set_PulseSpel_var('h',100)
        #     hw_log.warning('NOT USING PULSE_SPEL FOR FIELDSWEEP')
        #     self.api.set_Acquisition_mode(1) 
        #     self.api.run_exp()
        #     for i in range (120):
        #         if self.isrunning():
        #             time.sleep(1)
        #         else:
        #             break
        #     self.api.set_Acquisition_mode(0)
        #     self.api.set_xaxis_quanity('Magnetic Field')
            self.api.run_exp()
        pass

    
    def tune_rectpulse(self,*,tp, freq, B, reptime, shots=400):
        """Generates a rectangular pi and pi/2 pulse of the given length at 
        the given field position. This value is stored in the pulse cache. 

        Parameters
        ----------
        tp : float
            Pulse length in ns
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

        overflow_flag= True

        while overflow_flag:
            vg = self.api.get_video_gain()
            self.launch(amp_tune, "autoDEER_amptune")

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            dataset = dataset.epr.correctphase

            data = np.abs(dataset.data)
            scale_amp = np.around(dataset.pulse0_scale[data.argmax()].data,2)
            if scale_amp > 0.95:
                raise RuntimeError("Not enough power avaliable.")
            
            if scale_amp == 0:
                warnings.warn("Pulse tuned with a scale of zero!")
            p90 = amp_tune.pulses[0].copy(
                scale=scale_amp)
            
            p180 = amp_tune.pulses[1].copy(
                scale=scale_amp)
            
            # Now tune the corresponding videoGain

            amp_tune_short = amp_tune.copy()
            scale = Parameter("scale",scale_amp,dim=4,step=0.01)
            amp_tune_short.pulses[0].scale = scale
            amp_tune_short.pulses[1].scale = scale
            amp_tune_short.evolution([scale])

            self.launch(amp_tune_short, "autoDEER_amptune")

            time.sleep(3)
            self.terminate()
            time.sleep(2)
            while self.api.is_exp_running():
                time.sleep(1)

            self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
            TimeBase = self.api.hidden['specJet.TimeBase'].value 
            self.api.hidden['specJet.NoOfPoints'].value = int(512/TimeBase)
            self.api.hidden['specJet.NoOfAverages'].value = 1
            self.api.hidden['specJet.NOnBoardAvgs'].value = 1
            if not self.api.hidden['specJet.AverageStart'].value:
                self.api.hidden['specJet.AverageStart'].value = True
            time.sleep(2)
            limit = np.abs(self.api.hidden['specjet.DataRange'][0])
            data = get_specjet_data(self)
            
            if (np.abs(data.real).max() > 0.95*limit) or (np.abs(data.imag).max() > 0.95*limit):
                self.api.set_video_gain(vg - 12)
                overflow_flag= True
                print(f'overflow, videogain now {vg - 12}')
            elif (np.abs(data.real).max() < 0.02*limit) and (np.abs(data.imag).max() < 0.02*limit):
                self.api.set_video_gain(vg + 6)
                overflow_flag= True
                print(f'underflow, videogain now {vg + 6}')
            else:
                overflow_flag = False

        return p90, p180
    
    def tune_pulse(self, pulse, mode, freq, B , reptime, shots=400,tp=12):
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
        tp: float
            The length of the pi/2 pulse used for Hahn Echo, by default 12 ns. The pi pulse will be twice the length

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

        pi2_pulse, pi_pulse = self.tune_rectpulse(tp=tp, B=B, freq=c_frq, reptime=reptime)
        
        if mode == "amp_hahn":
            amp_tune =HahnEchoSequence(
                B=B, freq=freq, 
                reptime=reptime, averages=1, shots=shots,
                pi2_pulse = pulse, pi_pulse=pi_pulse
            )

            scale = Parameter('scale',0,unit=None,step=0.02, dim=51, description='The amplitude of the pulse 0-1')
            amp_tune.pulses[0].scale = scale

            amp_tune.evolution([scale])

            self.launch(amp_tune, "autoDEER_amptune")

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            new_amp = np.around(dataset.pulse0_scale[dataset.data.argmax()].data,2)
            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
            return pulse

        elif mode == "amp_nut":

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

            scale = Parameter('scale',0,unit=None,step=0.02, dim=51, description='The amplitude of the pulse 0-1')
            nut_tune.pulses[0].scale = scale
            nut_tune.evolution([scale])

            self.launch(nut_tune, "autoDEER_amptune")

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            dataset = dataset.epr.correctphase
            data = dataset.data
            axis = scale.get_axis()
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

    def phasetune_pulse(self, pulse):
        # Bruker SpecJet-I has a phase issue. The pulse is observed and optimised through the transmision mode.

        if isinstance(pulse, Detection):
            raise RuntimeError("Detection pulses cannot be phase tuned.")
        
        current_B = self.api.get_field()
        current_freq = self.api.get_counterfreq()

        test_seq = Sequence(name='test_seq',B=current_B,freq=current_freq,reptime=250,averages=1,shots=20)
        test_seq.addPulse(Detection(tp=1024,t=0))
        test_seq.addPulse(pulse.copy(t=200))
        test_seq.evolution([])

        self.launch(test_seq, "autoDEER_phasetune")
        time.sleep(1)
        self.terminate()

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
        self.api.hidden['specJet.NoOfAverages'].value = 20
        self.api.hidden['specJet.NOnBoardAvgs'].value = 20
        self.api.hidden['specJet.NoOfPoints'].value = 1024


        



        
    def isrunning(self) -> bool:
        if self.tuning:
            return True
        if self.bg_thread is None:
            return self.api.is_exp_running()
        else:
            return self.bg_thread.running()

    def terminate(self,now=False) -> None:
        self.tuning = False
        if self.bg_thread is None:
            if now:
                return self.api.abort_exp()
            else:
                return self.api.stop_exp()
        else:
            attempt = self.bg_thread.cancel()
            if not attempt:
                raise RuntimeError("Thread failed to be canceled!")
    

    # Override the terminate method to allow for a more complex termination, that reduces the number of call to the XeprAPI
    def terminate_at(self, criterion, test_interval=2, keep_running=True, verbosity=0,autosave=True):
        """Terminates the experiment upon a specific condition being
        satisified. 

        Parameters
        ----------
        criterion : _type_
            The criteria to be tested.
        test_interval : int, optional
            How often should the criteria be tested in minutes, by default 10.
        keep_running : bool, optional
            If True, an error will not be raised if the experiment finishes before the criteria is met, by default True.
        verbosity : int, optional
            The verbosity level, by default 0. 
        autosave : bool, optional
            If True, the data will be autosaved, by default True.

        """



        test_interval_seconds = test_interval * 60
        condition = False
        last_scan = 0

        seq_time_scan = self.cur_exp._estimate_time() / self.cur_exp.averages.value


        while not condition:
            
            time.sleep(np.max[seq_time_scan/4,10])

            if not self.isrunning():
                if keep_running:
                    self.terminate()
                    return None
                else:
                    msg = "Experiments has finished before criteria met."
                    raise RuntimeError(msg)

            start_time = time.time()

            data = self.acquire_dataset(after_scan=last_scan+1)
            if autosave:
                self.log.debug(f"Autosaving to {os.path.join(self.savefolder,self.savename)}")
                data.to_netcdf(os.path.join(self.savefolder,self.savename),engine='h5netcdf',invalid_netcdf=True)

            try:
                # nAvgs = data.num_scans.value
                nAvgs = data.attrs['nAvgs']

            except AttributeError or KeyError:
                self.log.warning("WARNING: Dataset missing number of averages(nAvgs)!")
                nAvgs = 1
            finally:
                if nAvgs < 1:
                    time.sleep(seq_time_scan/2)  # TODO: Replace with single scan time
                    continue
                elif nAvgs <= last_scan:
                    time.sleep(seq_time_scan/2)
                    continue    
            last_scan = nAvgs
            if verbosity > 0:
                print("Testing")

            if isinstance(criterion,list):
                conditions = [crit.test(data, verbosity) for crit in criterion]
                condition = any(conditions)

            else:
                condition = criterion.test(data, verbosity)

            if not condition:
                end_time = time.time()
                if (end_time - start_time) < test_interval_seconds:
                    if verbosity > 0:
                        print("Sleeping")
                    time.sleep(test_interval_seconds - (end_time - start_time))
        
        if isinstance(criterion,list):
            for i,crit in enumerate(criterion):
                if conditions[i]:
                    if callable(crit.end_signal):
                        crit.end_signal()
                
        else:
            if callable(criterion.end_signal):
                criterion.end_signal()
        
        
        self.terminate()
        pass
            
    def calc_d0(self):
        """
        This creates an initial guess for d0. A better estimate can only be found after the field sweep. 
        """
        hw_log.info('Calcuating d0')
        hw_log.debug('Setting Detection = TM')
        self.api.hidden['Detection'].value = 'TM'
        B = self.api.get_field()
        freq = self.api.get_counterfreq()

        self.api.set_attenuator('+<x>',100)

        d0=0
        self.d0=d0

        seq = Sequence(name='single_pulse',B=B,freq=freq,reptime=3e3,averages=1,shots=20)
        det_tp = Parameter('tp',value=16,dim=4,step=0)
        seq.addPulse(RectPulse(tp=det_tp,t=0,flipangle=np.pi,scale=1))
        seq.addPulse(Detection(tp=16,t=d0))

        seq.evolution([det_tp])
        self.launch(seq,savename='test',tune=False)
        self.terminate(now=True)
        time.sleep(3)

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
        self.api.hidden['specJet.NoOfAverages'].value = 20
        self.api.hidden['specJet.NOnBoardAvgs'].value = 20
                
        if not self.api.hidden['specJet.AverageStart'].value:
            self.api.hidden['specJet.AverageStart'].value = True

        # Set the time base of specjet to a minimum
        time_base = self.api.hidden['specJet.TimeBase'].aqGetParMinValue()
        self.api.hidden['specJet.TimeBase'].value = time_base
        n_points = int(1024/time_base)
        self.api.hidden['specJet.NoOfPoints'].value = n_points
        time.sleep(3)
        
        
        optimal = False
        while not optimal:
            max_value = np.abs(get_specjet_data(self)).max()
            y_max = np.abs(self.api.hidden['specjet.DataRange'][0])
            vg  =self.api.get_video_gain()
            vg_step = self.api.get_video_gain_step()
            if max_value > 0.7* y_max:
                self.api.set_video_gain(vg - vg_step)
                time.sleep(0.5)
            elif max_value < 0.3* y_max:
                self.api.set_video_gain(vg + vg_step)
                time.sleep(0.5)
            else:
                optimal=True
        
        specjet_data = np.abs(get_specjet_data(self))

        calc_d0 = d0  + self.api.hidden['specJet.Abs1Data'][specjet_data.argmax()]

        d0 = calc_d0 - 256

        seq = Sequence(name='single_pulse',B=B,freq=freq,reptime=3e3,averages=1,shots=20)
        det_tp = Parameter('tp',value=16,dim=4,step=0)
        seq.addPulse(RectPulse(tp=det_tp,t=0,flipangle=np.pi,scale=1))
        seq.addPulse(Detection(tp=16,t=d0))

        seq.evolution([det_tp])
        self.launch(seq,savename='test',tune=False)
        self.terminate()

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
                
        if not self.api.hidden['specJet.AverageStart'].value:
            self.api.hidden['specJet.AverageStart'].value = True

        self.api.hidden['specJet.NoOfPoints'].value = int(512/time_base)
        time.sleep(3)
        specjet_data = np.abs(get_specjet_data(self))
        
        calc_d0 = d0 + self.api.hidden['specJet.Abs1Data'][specjet_data.argmax()]

        self.d0 = calc_d0 + 64 # 64ns added to compensate for hahn echo center in field sweep
        
        hw_log.info(f"d0 set to {self.d0}")
        self.api.hidden['Detection'].value = 'Signal'

    def calc_d0_from_Hahn_Echo(self, B=None, freq=None):
        
        B = self.api.get_field()
        freq = self.api.get_counterfreq()

        if B is not None:
            self.api.set_field(B)
            time.sleep(5)
        else:
            B = self.api.get_field()
        if freq is not None:
            self.api.set_freq(freq)
        else:
            freq = self.api.get_counterfreq()
        
        d0 = self.d0
        # self.api.set_PulseSpel_var('d0',d0)
        # self.api.run_exp()
        # self.api.abort_exp()
        while self.api.is_exp_running():
            time.sleep(1)

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
        self.api.hidden['specJet.NoOfAverages'].value = 20
        self.api.hidden['specJet.NOnBoardAvgs'].value = 20
                
        if not self.api.hidden['specJet.AverageStart'].value:
            self.api.hidden['specJet.AverageStart'].value = True

        time_base = self.api.hidden['specJet.TimeBase'].aqGetParMinValue()
        self.api.hidden['specJet.TimeBase'].value = time_base

        self.api.hidden['specJet.NoOfPoints'].value = int(512/time_base)

        optimal = False
        while not optimal:
            max_value = np.abs(get_specjet_data(self)).max()
            y_max = np.abs(self.api.hidden['specjet.DataRange'][0])
            vg  =self.api.get_video_gain()
            vg_step = self.api.get_video_gain_step()
            if max_value > 0.7* y_max:
                self.api.set_video_gain(vg - vg_step)
                time.sleep(0.5)
            elif max_value < 0.3* y_max:
                self.api.set_video_gain(vg + vg_step)
                time.sleep(0.5)
            else:
                optimal=True
        
        specjet_data = np.abs(get_specjet_data(self))
        calc_d0 = d0 - 64 + self.api.hidden['specJet.Abs1Data'][specjet_data.argmax()]

        self.d0 = calc_d0 
        hw_log.info(f"d0 set to {self.d0}")
        return self.d0

# =============================================================================


def test_if_MPFU_compatability(seq):
    table = seq.progTable
    if 'freq' in table['Variable']:
        return False
    elif np.unique(table['axID']).size > 2:
        return False 
    else:
        return True