import numpy as np
import os
import copy
import time
import datetime
from pathlib import Path
import logging
import yaml
import warnings

from pyepr.utils import round_step
from pyepr.config import get_waveform_precision
import pyepr.pulses as pulses
from pyepr.classes import Parameter
from pyepr.sequences import HahnEchoSequence



# =============================================================================


class Interface:
    """Represents the interface connection from autoEPR to the spectrometer.
    """

    def __init__(self,config_file:dict=None,log=None) -> None:
        """
        Parameters
        ----------
        config_file : dict or str or Path, optional
            The configuration file or dict for the spectrometer interface, by default None. If None, a default configuration will be used.
        log : logging.Logger, optional  
            The logger to be used, by default None. If None, a default logger will be created.
        """
        if isinstance(config_file, (str,Path)):
            with open(config_file, 'r') as f:
                config_file = yaml.safe_load(f)
        
        self.config = config_file if isinstance(config_file, dict) else {"Spectrometer":{"Bridge":{}}}
        
        self.pulses = {}
        self.savefolder = str(Path.home())
        self.savename = ""
        if log is None:
            self.log = logging.getLogger('interface')
        else:
            self.log = log
        self.resonator = None
        if self.config != {}:
            self.amp_nonlinearity = self.config["Spectrometer"]["Bridge"].get('Amplifier Non-Linearity',None)
        else:
            self.amp_nonlinearity = None

        self.AWG=True
        pass

    def connect(self) -> None:
        pass

    def acquire_dataset(self, data):
        """
        Acquires the dataset.
        """

        # data.sequence = self.cur_exp
        data.attrs['time'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        return data
    

    def launch(self, sequence, savename: str):
        """Launches the experiment and initialises autosaving.

        Parameters
        ----------
        sequence : Sequence
            The sequence to be launched
        savename : str
            The savename for this measurement. A timestamp will be added to the value.
        """
        timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
        self.savename=timestamp + savename + '.h5'
        pass

    def isrunning(self) -> bool:
        return False
    
    def rescale(self,scale: float) -> float:
        """Rescales the scale factor to account for applifier non-linearity.

        Parameters
        ----------
        scale : float
            The scale factor to be rescaled.

        Returns
        -------
        float
            The rescaled value.
        """

        if self.amp_nonlinearity is None:
            print("WARNING: No amplifier non-linearity defined. Using linear scaling.")
            return scale
        if isinstance(self.amp_nonlinearity, list): # Assume polyniomial coefficients 
            coeff = copy.copy(self.amp_nonlinearity)
            coeff[-1] -= scale  # Set the last coefficient to the negative of the scale factor
            roots = np.roots(coeff) # Check if the coefficients are valid
            real_roots = roots[np.isreal(roots)].real

            if np.all(real_roots < 0):
                new_scale = scale
                print("Warning: All roots are negative, setting scale to orginal")
            elif len(real_roots[(real_roots >= 0) & (real_roots <= 1)]) == 0:
                new_scale = 1.0
            else:
                valid = real_roots[(real_roots >= 0) & (real_roots <= 1)]
                new_scale = valid[0]

            new_scale = np.clip(new_scale, 0, 1)
            return new_scale

    def terminate(self) -> None:
        """
        Terminates the experiment immediately. 
        """
        pass

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


        while not condition:
            
            time.sleep(10) # TODO: Replace with half sequence time 

            if not self.isrunning():
                if keep_running:
                    self.terminate()
                    return None
                else:
                    msg = "Experiments has finished before criteria met."
                    raise RuntimeError(msg)

            start_time = time.time()
            data = self.acquire_dataset()
            if autosave:
                self.log.debug(f"Autosaving to {os.path.join(self.savefolder,self.savename)}")
                # data.to_netcdf(os.path.join(self.savefolder,self.savename),engine='h5netcdf',invalid_netcdf=True)
                data.epr.save(os.path.join(self.savefolder,self.savename))

            try:
                # nAvgs = data.num_scans.value
                nAvgs = data.attrs['nAvgs']

            except AttributeError or KeyError:
                self.log.warning("WARNING: Dataset missing number of averages(nAvgs)!")
                nAvgs = 1
            finally:
                if nAvgs < 1:
                    time.sleep(30)  # TODO: Replace with single scan time
                    continue
                elif nAvgs <= last_scan:
                    time.sleep(30)
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

    def tune_RectPulse(self,freq, B, reptime, tp, offset_freq=0, shots=200, verbosity=0,**kwargs):
        """
        Tunes a pair of equal power pi/2, pi rectangular pulses though an amplitude sweep. 
        
        Parameters
        ----------
        freq : float
            The central frequency of the experiment in GHz.
        B : float
            The magnetic field in Gauss.
        reptime : float
            The repetition time in us.
        tp : float
            The length of the pi/2 pulse in ns. The pi pulse will be half this length. 
            This will be rounded to the nearest 2*waveform_precision to ensure that the pi pulse is of the correct length.
        offset_freq : float, optional
            The offset frequency of the rectangular pulse in GHz. This is used to shift the pulse away from the central frequency of the experiment. Default is 0.
        shots : int, optional
            The number of shots to be used for the frequency sweep, by default 200.
        verbosity : int, optional
            The verbosity level, by default 0. If > 0, the scale values will be printed.
        **kwargs : dict
            - step : float, optional
                The step size for the amplitude parameter, by default 0.02.
            - dim : int, optional
                The dimension of the amplitude parameter, by default 50.
            - tau : float, optional
                The delay between the pi/2 and pi pulses in ns, by default 500.
        ** Interface specific launch commands

        Returns
        -------
        p90 : pulses.RectPulse
            The tuned pi/2 pulse.
        p180 : pulses.RectPulse
            The tuned pi pulse.
        """
        step = kwargs.pop('step', 0.02)
        dim = kwargs.pop('dim', 50)
        start = kwargs.pop('start', 0.0)
        tau = kwargs.pop('tau', 500)

        amp = Parameter(name="scale", value=start, step=step, dim=dim)
        p90 = pulses.RectPulse(tp=tp,freq=offset_freq,flipangle=np.pi/2,scale=amp)
        p180 = pulses.RectPulse(tp=tp*2,freq=offset_freq,flipangle=np.pi,scale=amp)
        RectEchoTuneSeq =HahnEchoSequence(
                    B=B, freq=freq, reptime=reptime, averages=1, shots=shots, tau=tau,
                    pi2_pulse=p90, pi_pulse=p180,
                    name='RectPulse Echo Tune Sequence')

        RectEchoTuneSeq.evolution([amp])

        self.launch(RectEchoTuneSeq, 
                    savename=f"RectEchoTune",
                    **kwargs)
        time.sleep(5)
        while self.isrunning():
            time.sleep(1)

        dataset = self.acquire_dataset(filter_type='cheby',filter_width=p90.bandwidth.value*1e3)

        if np.any(dataset.data == 0):
            warnings.warn("Zero values found in dataset. This may indicate an error in acquisition.")
            time.sleep(5)
            while self.isrunning():
                time.sleep(2)
            dataset = self.acquire_dataset(downconvert=True,reduce=True,filter_type='boxcar',filter_width=250)

        dataset = dataset.epr.correctphase
        data = np.abs(dataset.data)
        scale = np.around(dataset.pulse0_scale[data.argmax()].data,2)
        if verbosity > 0:
            print(f"Scale Value: {scale}")
        if scale > 0.95:
            raise RuntimeError("Not enough power avaliable.")

        if scale == 0:
            warnings.warn("Pulse tuned with a scale of zero!")
            print("Pulse tuned with a scale of zero!")
        p90 = p90.copy(scale=scale)
        p180 = p180.copy(scale=scale)

        return p90, p180



    def tune_GaussianPulse(self,freq, B, reptime, tp, offset_freq=0, shots=200, verbosity=0,**kwargs):
        """
        Tunes a pair of equal power pi/2, pi Gaussian pulses though an amplitude sweep. 
        
        Parameters
        ----------
        freq : float
            The central frequency of the experiment in GHz.
        B : float
            The magnetic field in Gauss.
        reptime : float
            The repetition time in us.
        tp : float
            The length of the pi/2 pulse in ns. The pi pulse will be half this length. 
            This will be rounded to the nearest 2*waveform_precision to ensure that the pi pulse is of the correct length.
        offset_freq : float, optional
            The offset frequency of the Gaussian pulse in GHz. This is used to shift the pulse away from the central frequency of the experiment. Default is 0.
        shots : int, optional
            The number of shots to be used for the frequency sweep, by default 200.
        verbosity : int, optional
            The verbosity level, by default 0. If > 0, the scale values will be printed.
        **kwargs : dict
            - step : float, optional
                The step size for the amplitude parameter, by default 0.02.
            - dim : int, optional
                The dimension of the amplitude parameter, by default 50.
            - tau : float, optional
                The delay between the pi/2 and pi pulses in ns, by default 500.
        ** Interface specific launch commands

        Returns
        -------
        p90 : pulses.GaussianPulse
            The tuned pi/2 pulse.
        p180 : pulses.GaussianPulse
            The tuned pi pulse.

        """

        if not self.AWG:
            raise RuntimeError("Only AWG based spectrometers can tune Gaussian pulses.")

        step = kwargs.pop('step', 0.02)
        dim = kwargs.pop('dim', 50)
        start = kwargs.pop('start', 0.0)
        tau = kwargs.pop('tau', 500)

        amp = Parameter(name="scale", value=start, step=step, dim=dim)
        p90 = pulses.GaussianPulse(tp=tp,freq=offset_freq,flipangle=np.pi/2,scale=amp)
        p180 = pulses.GaussianPulse(tp=tp*2,freq=offset_freq,flipangle=np.pi,scale=amp)
        GaussianEchoTuneSeq =HahnEchoSequence(
                    B=B, freq=freq, reptime=reptime, averages=1, shots=shots, tau=tau,
                    pi2_pulse=p90, pi_pulse=p180,
                    name='GaussianPulse Echo Tune Sequence')

        GaussianEchoTuneSeq.evolution([amp])

        self.launch(GaussianEchoTuneSeq, 
                    savename=f"GaussEchoTune",
                    **kwargs)
        time.sleep(5)
        while self.isrunning():
            time.sleep(1)

        dataset = self.acquire_dataset(filter_type='cheby',filter_width=p90.bandwidth.value*1e3)

        if np.any(dataset.data == 0):
            warnings.warn("Zero values found in dataset. This may indicate an error in acquisition.")
            time.sleep(5)
            while self.isrunning():
                time.sleep(2)
            dataset = self.acquire_dataset(downconvert=True,reduce=True,filter_type='boxcar',filter_width=250)

        dataset = dataset.epr.correctphase
        data = np.abs(dataset.data)
        scale = np.around(dataset.pulse0_scale[data.argmax()].data,2)
        if verbosity > 0:
            print(f"Scale Value: {scale}")
        if scale > 0.95:
            raise RuntimeError("Not enough power avaliable.")

        if scale == 0:
            warnings.warn("Pulse tuned with a scale of zero!")
            print("Pulse tuned with a scale of zero!")
        p90 = p90.copy(scale=scale)
        p180 = p180.copy(scale=scale)

        return p90, p180

    
    def tune_KBBPulse(self,*, freq, B, reptime, tp, init_freq, final_freq, shots=200, verbosity=0,**kwargs):
        """
        Tunes a Kutz-Boehlen-Bodenhausen linear chirp echo. This is a three step optimisation:

        1. Both pulses amplitudes are swept with a 1:2 ratio.
        2. Pulse 0 is fixed using the optimal amplitude from step 1 and pulse 1 is swept.
        3. Pulse 1 is fixed using the optimal amplitude from step 2 and pulse 0 is swept.

        Parameters
        ----------
        freq : float
            The central frequency of the experiment in GHz.
        B : float
            The magnetic field in Gauss.
        reptime : float
            The repetition time in us.
        tp : float
            The length of the pi/2 pulse in ns. The pi pulse will be half this length. 
            This will be rounded to the nearest 2*waveform_precision to ensure that the pi pulse is of the correct length.
        init_freq : float
            The initial frequency of the frequency sweep in GHz.
        final_freq : float
            The final frequency of the frequency sweep in GHz.
        shots : int, optional
            The number of shots to be used for the frequency sweep, by default 200.
        verbosity : int, optional
            The verbosity level, by default 0. If > 0, the scale values will be printed.
        **kwargs : dict
            - step : float, optional
                The step size for the amplitude parameter, by default 0.01.
            - dim : int, optional
                The dimension of the amplitude parameter, by default 50.
            - tau : float, optional
                The delay between the pi/2 and pi pulses in ns, by default 500.
        ** Interface specific launch commands

        Returns
        -------
        p90 : pulses.ChirpPulse
            The tuned pi/2 pulse.
        p180 : pulses.ChirpPulse
            The tuned pi pulse.
        """
        if not self.AWG:
            raise RuntimeError("Only AWG based spectrometers can tune KBB pulses.")

        waveform_precision = get_waveform_precision()
        tp = round_step(tp,waveform_precision*2) # Round to nearest 2*waveform_precision to ensure that the pi pulse is an integer multiple of waveform_precision

        step = kwargs.pop('step', 0.01)
        dim = kwargs.pop('dim', 50)
        start = kwargs.pop('start', 0.0)
        amp = Parameter(name="Scale", value=start, step=step, dim=dim)

        p90 = pulses.ChirpPulse(
            tp=tp,init_freq=init_freq,final_freq=final_freq,flipangle=np.pi/2,scale=amp)
        p180 = pulses.ChirpPulse(
            tp=tp*0.5,init_freq=init_freq,final_freq=final_freq,flipangle=np.pi,scale=amp)


        tau = kwargs.pop('tau', 500)

        ChirpEchoTuneSeq = HahnEchoSequence(
            freq=freq,B=B, reptime=reptime, shots=shots, averages=1, pi2_pulse=p90, 
            pi_pulse=p180, tau=tau)
        ChirpEchoTuneSeq.evolution([amp])


        self.launch(ChirpEchoTuneSeq, 
                    savename=f"chirpechotune_step1",
                    **kwargs)
        time.sleep(5)
        while self.isrunning():
            time.sleep(1)

        dataset = self.acquire_dataset(filter_type='cheby',filter_width=p90.bandwidth.value*1e3)

        scale = np.around(dataset.pulse0_scale[np.abs(dataset.epr.correctphase).argmax()].data,2)

        if verbosity > 0:
            print(f"Initial Scale Value: {scale}")

        ChirpEchoTuneSeq.pulses[0].scale = Parameter('scale',value=scale)
        ChirpEchoTuneSeq.pulses[1].scale = amp*2
        ChirpEchoTuneSeq.evolution([amp])

        self.launch(ChirpEchoTuneSeq, 
                            savename=f"chirpechotune_step2",
                            **kwargs)
        
        while self.isrunning():
            time.sleep(1)
        
        dataset = self.acquire_dataset(filter_type='cheby',filter_width=p90.bandwidth.value*1e3)
        
        scale1 = np.around(dataset.pulse1_scale[np.abs(dataset.epr.correctphase).argmax()].data,2)

        if verbosity > 0:
            print(f"Improved pulse1 scale: {scale1}")

        ChirpEchoTuneSeq.pulses[0].scale = amp
        ChirpEchoTuneSeq.pulses[1].scale = Parameter('scale',value=scale)
        ChirpEchoTuneSeq.evolution([amp])

        self.launch(ChirpEchoTuneSeq, 
                                    savename=f"chirpechotune_step3",
                                    **kwargs)

        while self.isrunning():
                    time.sleep(1)
                
        dataset = self.acquire_dataset(filter_type='cheby',filter_width=p90.bandwidth.value*1e3)
        
        scale0 = np.around(dataset.pulse0_scale[np.abs(dataset.epr.correctphase).argmax()].data,2)

        if verbosity > 0:
            print(f"Improved pulse0 scale: {scale0}")

        p90 = pulses.ChirpPulse(
            tp=tp,init_freq=init_freq,final_freq=final_freq,flipangle=np.pi/2,scale=scale0)
        p180 = pulses.ChirpPulse(
            tp=tp*0.5,init_freq=init_freq,final_freq=final_freq,flipangle=np.pi,scale=scale1)

        return p90, p180

