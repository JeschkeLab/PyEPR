from pyepr.classes import Parameter
from pyepr.pulses import Pulse, Detection, Delay, RectPulse
from pyepr.utils import val_in_ns, val_in_us, add_phaseshift, _gen_ESEEM
from pyepr.fieldsweep_analysis import create_Nmodel
import pyepr.pulses as ad_pulses
import numpy as np
import matplotlib.pyplot as plt
import json
from itertools import product
import copy
from pyepr.utils import build_table, autoEPRDecoder
import json
from pyepr import __version__
import uuid
import base64
import numbers
import warnings

class Sequence:
    """
    Represents an experimental pulse sequence.
    """

    def __init__(
            self, *, name, B, reptime, averages, shots,LO=None, freq=None, **kwargs) -> None:
        """Represents an experimental pulse sequence.

        Parameters
        ----------
        name : str
            The name of this pulse sequence
        B : float
            The magnetic field for this sequence in Gauss.
        freq : float
            The central frequency of this sequence. I.e. The frequnecy at which
            a zero offset pulse is at. 
        reptime : float
            The shot repetition time in us.
        averages : int
            The number of scans to be accumulated.
        shots : itn
            The number of shots per point.
        LO : float
            The now deprecated local oscillator frequency. freq should be used.
        """

        self.pulses = []
        self.num_pulses = len(self.pulses)
        self.axes_uuid = []
        self.reduce_uuid = []



        if isinstance(B, Parameter):
            self.B = B.copy()
        else:
            self.B = Parameter(
                "B", B, "Gauss",
                "The static B0 field for the experiment")
        
        if freq is not None:
            self.freq = Parameter(
                "freq", freq, "GHz",
                "The central frequency of the sequence")
        elif LO is not None:
            self.freq = Parameter(
                "freq", LO, "GHz",
                "The central frequency of the sequence")
            warnings.warn("LO is deprecated, please use freq instead", DeprecationWarning, stacklevel=2)
        else:
            raise ValueError("Either freq or LO must be specified")
        
        if isinstance(reptime, Parameter):
            self.reptime = reptime.copy()
        else:
            self.reptime = Parameter(
                "reptime", reptime, "us",
                "The shot repetition time")
        
        self.averages = Parameter(
            "averages", averages, "None",
            "The number of averages to perform.")
        
        self.shots = Parameter(
            "shots", shots, "None",
            "The number of shots per scan.")

        if "det_window" in kwargs:
            self.det_window = Parameter(
                "det_window", kwargs["det_window"], "None",
                "The length of the default detection gate"
            )
        else:
            self.det_window = Parameter(
                "det_window", 128, "None",
                "The length of the default detection gate"
            )

            
        self.name = name
        self.progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": [], "uuid": [], "reduce": []}
        pass

    def plot(self) -> None:

        pass

    def plot_pulse_exc(self, FieldSweep=None, ResonatorProfile=None):

        if FieldSweep is not None:
            data = FieldSweep.data
            data /= np.max(np.abs(data))
            ax = FieldSweep.fs_x
            plt.plot(ax, data, label="Field Sweep", color='r')
            
        if ResonatorProfile is not None:
            data = ResonatorProfile.prof_data / ResonatorProfile.prof_data.max()
            plt.plot(ResonatorProfile.prof_frqs-self.LO.value, data, label="Resonator Profile", color='k')


        def check_array_new(list,array, return_pos=False):
            test = True
            for i, ele in enumerate(list):
                if np.allclose(
                    ele,array,rtol=1e-03, atol=1e-03, equal_nan=True):
                    test = False
                    if return_pos:
                        return i
            return test

        axs = []
        fts = []
        labels = []
        for i, pulse in enumerate(self.pulses):
            if type(pulse) == Delay:
                continue
            elif type(pulse) == Detection:
                continue
            ax, ft = pulse._calc_fft(1e4)

            if check_array_new(fts, ft):
                axs.append(ax)
                fts.append(ft)
                labels.append([i])
            else:
                j = check_array_new(fts, ft, return_pos=True)
                labels[j].append(i)

        hatches = ['/', '\\', '|', '-', 'o', '.']
        
        
        for i in range(0,len(axs)):
            ft = np.abs(fts[i])
            ft /= np.max(ft)
            plt.fill(axs[i],ft, 
                     label=f"Pulse id={str(labels[i])[1:-1]}",
                     hatch=hatches.pop(), alpha=0.3)
        
        plt.xlabel("Frequency (GHz)")
        plt.legend()

        # Set correct axis, find mim and max freq
        
        for i, pulse in enumerate(self.pulses):
            if hasattr(pulse, "freq"):
                p_min = pulse.freq.value
                p_max = pulse.freq.value
            elif hasattr(pulse, "init_freq") & hasattr(pulse, "final_freq"):
                p_min = pulse.init_freq.value
                p_max = pulse.final_freq.value
            elif hasattr(pulse, "init_freq") & hasattr(pulse, "BW"):
                p_min = pulse.init_freq.value
                p_max = p_min + pulse.BW.value
            elif hasattr(pulse, "final_freq") & hasattr(pulse, "BW"):
                p_max = pulse.init_freq.value
                p_min = p_max - pulse.BW.value
            # Check p_min < p_max
            if p_min > p_max:
                p_temp = p_min
                p_min = p_max
                p_max = p_temp
            if i == 0:
                min_frq = p_min
                max_frq = p_max
            else:
                if min_frq > p_min:
                    min_frq = p_min
                if max_frq < p_max:
                    max_frq = p_max
        
        min_frq = np.floor(min_frq*100)/100 - 0.2
        max_frq = np.ceil(max_frq*100)/100 + 0.2
        plt.xlim(min_frq, max_frq)
        
        
        pass

    def addPulse(self, pulse):
        """Adds a pulse to the next position in the sequence.

        Parameters
        ----------
        pulse : Pulse
            The object describing the pulse.
        """
        if type(pulse) == Pulse or issubclass(type(pulse), Pulse):
            self.pulses.append(pulse)
        
        elif type(pulse) == list:
            for el in pulse:
                self.pulses.append(el)
        self.num_pulses = len(self.pulses)
        self._buildPhaseCycle()
        self._estimate_time()

    def _estimate_time(self):
        """
        Calculates the estimated experiment time in seconds.
        """
        self._buildPhaseCycle()
        acqs = self.averages.value * self.shots.value * self.pcyc_dets.shape[0]
        if hasattr(self,'evo_params'):
            acqs *= np.prod([np.prod(param.dim) for param in self.evo_params])        
        time = acqs * self.reptime.value * 1e-6

        self.time = Parameter(name="time", value=f"{(time // 3600):.0f}:{(time % 3600) // 60:.0f}:{(time % 60):.0f}", unit="HH:MM:SS",
        description="Estimated sequence run time")
        return time

    def _buildPhaseCycle(self):
        # Identify pulses which are phase cycled

        pcyc_pulses = []
        pulse_cycles = []
        det_cycles = []
        for ix, pulse in enumerate(self.pulses):
            if pulse.pcyc is not None:
                pcyc_pulses.append(ix)
                # pulse_cycles.append(np.array(pulse.pcyc[0]))
                # det_cycles.append(np.array(pulse.pcyc[1]))
                pulse_cycles.append(pulse.pcyc["Phases"])
                det_cycles.append(pulse.pcyc["DetSigns"])

        self.pcyc_cycles = np.array(list(product(*pulse_cycles)))
        self.pcyc_dets = np.array(list(product(*det_cycles))).prod(axis=1)
        self.pcyc_vars = pcyc_pulses

        # # Build expanded phase cycle
        # func = lambda x: np.arange(0, len(x))
        # map(func, pulse_cycles)
        # n_pulses = len(pulse_cycles)

        # m = list(map(func, pulse_cycles))
        # grids = np.meshgrid(*m, indexing='ij')
        # expanded_cycles = []
        # expanded_dets = []

        # for i in range(0, n_pulses):
        #     expanded_cycles.append(
        #         pulse_cycles[i][grids[i].flatten(order='F')])
        #     expanded_dets.append(det_cycles[i][grids[i].flatten(order='F')])

        # self.pcyc_vars = pcyc_pulses
        # self.pcyc_cycles = np.stack(expanded_cycles)
        # self.pcyc_dets = np.prod(np.stack(expanded_dets), axis=0)
        # self.pcyc_n = self.pcyc_cycles.shape[1]
        
        return self.pcyc_vars, self.pcyc_cycles, self.pcyc_dets

    def evolution(self, params, reduce=[]):
        """
        Sets what parameters are being evolved in the sequence, and which are
        being automatically.

        `self.evo_params = params`

        Parameters
        ----------
        params : list
            A list of Parameter objects which are being evolved. Each every 
            entry in the list will be a new axis in the sequence. Only one parameter
            per axis should be specified.
        reduce : list
            A list of Parameter objects which are being reduced. These are the
            parameters which are being averaged over. These parameters should
            also be in the params list.
        
        Returns
        -------
        progTable : dict
            A dictionary containing the progression of the sequence.
        
        
        """
        self.evo_params = params
        self.axes_uuid = [param.uuid for param in params]
        self.reduce_uuid = [param.uuid for param in reduce]

        self._buildProgTable()


    def _buildProgTable(self):

        progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": [], "uuid": [], "reduce": []}
        
        for n, pulse in enumerate(self.pulses):
            table = pulse.build_table()
            for i in range(len(table["uuid"])):
                if table["uuid"][i] in self.axes_uuid:
                    progTable["axID"].append(self.axes_uuid.index(table["uuid"][i]))
                    progTable["uuid"].append(table["uuid"][i]) 
                    progTable["EventID"].append(n)
                    progTable["Variable"].append(table["Variable"][i])
                    progTable["axis"].append(table["axis"][i])
                    if table["uuid"][i] in self.reduce_uuid:
                        progTable["reduce"].append(True)
                    else:
                        progTable["reduce"].append(False)

        for var_name in vars(self):
            var = getattr(self, var_name)
            if type(var) is Parameter:
                if not var.is_static() and not var.virtual:
                    for i in range(len(var.axis)):
                        if var.axis[i]["uuid"] in self.axes_uuid:
                            progTable["axID"].append(self.axes_uuid.index(var.axis[i]["uuid"]))
                            progTable["EventID"].append(None)
                            progTable["Variable"].append(var_name) 
                            progTable["axis" ].append(var.axis[i]["axis"])
                            progTable["uuid"].append(var.axis[i]["uuid"]) 
                            if var.axis[i]["uuid"] in self.reduce_uuid:
                                progTable["reduce"].append(True)
                            else:
                                progTable["reduce"].append(False)
        self.progTable = progTable
        self._estimate_time()
        return self.progTable
    
    @property
    def seqtable_steps(self):
        if len(self.evo_params) > 0:
            return self.pcyc_dets.shape[0] * (len(self.pcyc_vars)+1) * np.prod([np.prod(param.dim) for param in self.evo_params])
    
    @property
    def shape(self):
        """
        Gives the shape of the sequence,excluding the number of shots. 
        Return
        ------
        list:
            [nAvgs, axes_dims, nPcyc]

        """
        nAvgs = self.averages.value
        nPcyc = self.pcyc_dets.shape[0]
        nAxes = len(self.evo_params)
        if nAxes > 0:
            axes_dim  = []
            for i in range(nAxes):
                axes_dim.append(self.evo_params[i].dim[0])
        else:
            axes_dim = [1]

        return [nAvgs]+  axes_dim +[nPcyc]



    def adjust_step(self,waveform_precision):
        """
        Adjust the step size of all axes and pulses to be an integer multiple of the waveform precision
        This is to ensure that the waveform is generated correctly by the specific AWG

        Parameters
        ----------
        waveform_precision : float
            The precision of the waveform in ns

        """

        for param in self.evo_params:
            param.adjust_step(waveform_precision)
        for pulse in self.pulses:
            pulse.t.adjust_step(waveform_precision)
            pulse.tp.adjust_step(waveform_precision)

        self._buildProgTable()
        return self

    def shift_detfreq_to_zero(self):
        det_pulse = None
        for pulse in self.pulses:
            if isinstance(pulse,Detection):
                det_pulse = pulse
        
        det_freq = det_pulse.freq.value
        self.freq.value -= det_freq
        for pulse in self.pulses:
            if hasattr(pulse,'freq'):
                pulse.freq.value -= det_freq
            if hasattr(pulse,'init_freq'):
                pulse.init_freq.value -= det_freq
            if hasattr(pulse,'final_freq'):
                pulse.final_freq.value -= det_freq
        return self
    
    def _checkRect(self) -> bool:
        """Checks if all the pulses in the sequence are rectangular.
        """
        test = True

        for pulse in self.pulses:
            if type(pulse) is not RectPulse:
                test = False

        return test

    def __str__(self):

        header = "#" * 79 + "\n" + "PyEPR Sequence Definition" + \
                 "\n" + "#" * 79 + "\n"

        # Sequence Parameters
        seq_param_string = "Sequence Parameters: \n"
        seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
            'Name', 'Value', 'Unit', 'Description')

        for param_key in vars(self):
            param = getattr(self, param_key)
            if type(param) is Parameter:
                if param.unit is None:
                    unit = ""
                else:
                    unit = param.unit
                if type(param.value) is str:
                    seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
                        param.name, param.value, unit, param.description)
                else:
                    seq_param_string += "{:<10} {:<12.5g} {:<10} {:<30} \n".format(
                        param.name, param.value, unit, param.description)
        
        # Pulses
        pulses_string = "\nEvents (Pulses, Delays, etc...): \n"

        params = ['iD', 't', 'tp', 'scale', 'type', 'Phase Cycle']
        params_widths = ["4", "8", "8", "8", "14", "40"]

        pulses_string += build_table(self.pulses, params, params_widths)

        def print_event_id(i):
            if prog_table["EventID"][i] is None:
                return "Seq"
            else:
                return str(prog_table["EventID"][i])

        def get_unit(i):
            pulse_num = prog_table["EventID"][i]
            if pulse_num is None:
                param = getattr(self, prog_table["Variable"][i])
            else:
                param = getattr(
                    self.pulses[pulse_num], prog_table["Variable"][i])
            
            if param.unit is None:
                return "None"
            else:
                return param.unit


        def test_unique_step(array):
            diffs = np.diff(array)-np.diff(array)[0]
            return np.isclose(diffs,np.zeros(diffs.shape)).all()

        # Progressive elements
        prog_string = "\nProgression: \n"
        if len(self.progTable["axID"]) >= 1:
            prog_string += "{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} \n".format(
                'Pulse', 'Prog. Axis', 'Parameter', 'Step', 'Dim', 'Unit')
        for i in range(0, len(self.progTable["axID"])):
            prog_table = self.progTable
            axis = prog_table["axis"][i]
            if test_unique_step(axis):
                step = np.unique(np.diff(axis))[0]
                fstring = "{:<10} {:<10} {:<10} {:<10.5g} {:<10} {:<10} \n"
            else:
                step = "Var"
                fstring = "{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} \n"
             
            prog_string += fstring.format(
                print_event_id(i), prog_table["axID"][i],
                prog_table["Variable"][i], step,
                prog_table["axis"][i].shape[0], get_unit(i))

        footer = "#" * 79 + "\n" +\
            f"Built by PyEPR Version: {__version__}" + "\n" + "#" * 79

        return header + seq_param_string + pulses_string + prog_string + footer

    def copy(self):
        return copy.deepcopy(self)
    
    def _to_dict(self):
        to_return = {"version": __version__, "type": "Sequence", "subclass": str(type(self))}

        for key, var in vars(self).items():
            if isinstance(var, Parameter):
                to_return[key] = var._to_dict()
            if key == "pulses":
                new_list = []
                for pulse in var:
                    new_list.append(pulse._to_dict())
                to_return[key] = new_list
            else:
                to_return[key] = var

        return to_return

    def _to_json(self):
        class autoEPREncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    if (len(obj) > 0 ) and isinstance(obj[0], str):
                        return list(obj)
                    data = np.ascontiguousarray(obj.data)
                    data_b64 = base64.b64encode(data)
                    return dict(__ndarray__=str(data_b64),
                                dtype=str(obj.dtype),
                                shape=obj.shape)
                if isinstance(obj, complex):
                    return str(obj)
                if isinstance(obj, numbers.Number):
                    return str(obj)
                if isinstance(obj, uuid.UUID):
                    return_dict = {"__uuid__": str(obj)}
                    return return_dict
                if isinstance(obj, Parameter):
                    return obj._to_dict()
                if isinstance(obj, Pulse):
                    return obj._to_dict()
                if isinstance(obj, Sequence):
                    return obj._to_dict()
                else:
                    return json.JSONEncoder.default(self, obj)
        
        return json.dumps(self._to_dict(), cls=autoEPREncoder, indent=4)

    def save(self, filename):
        """Save the sequence to a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the object cannot be serialized to JSON.

        Example
        -------
        >>> obj = Sequence()
        >>> obj.save("my_sequence.json")
        """
        
                
        with open(filename, "w") as f:
           f.write(self._to_json())

    @classmethod
    def _from_dict(cls, dct):
        name = dct["name"]
        B = Parameter._from_dict(dct["B"])
        freq = Parameter._from_dict(dct["freq"])
        reptime = Parameter._from_dict(dct["reptime"])
        averages = Parameter._from_dict(dct["averages"])
        shots = Parameter._from_dict(dct["shots"])
        new_sequence = cls(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages, shots=shots
        )
        for key, var in dct.items(): 
            if isinstance(var, dict) and ("type" in var):
                setattr(new_sequence, key, Parameter._from_dict(var))
            elif key == "pulses":
                for pulse in var:
                    if hasattr(ad_pulses,pulse['type']):
                        new_sequence.pulses.append(
                            getattr(ad_pulses,pulse['type'])._from_dict(pulse))
                    else:
                        new_sequence.pulses.append(
                            Pulse._from_dict(pulse))
            else:
                setattr(new_sequence, key, var)

        return new_sequence
    
    @classmethod
    def _from_json(cls, JSONstring):
        dct = json.loads(JSONstring, object_hook=autoEPRDecoder)
        return cls._from_dict(dct)
    
    @classmethod
    def load(cls, filename):
        """Load an object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        obj : Sequence
            The Sequence loaded from the JSON file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.

        Example
        -------
        >>> obj = Sequence.load("my_sequence.json")
        """
        with open(filename, "r") as f:
           file_buffer = f.read()
        return cls._from_json(file_buffer)

# =============================================================================
#                                Subclasses
# =============================================================================

class HahnEchoSequence(Sequence):
    """
    Represents a Hahn-Echo sequence. 
    """
    def __init__(self, *, B, freq, reptime, averages, shots, **kwargs) -> None:
        """Build a Hahn-Echo sequence using either rectangular pulses or
        specified pulses. By default no progression is added to this sequence.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        
        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        det_event : Pulse
            An autoEPR Pulse object describing the detection even. If
            not specified a standard detection event will be created instead,
            with an offset frequency of 0MHz. 

        """
        
        name = "HahnEchoSequence"
        
        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]
        if "det_event" in kwargs:
            self.det_event = kwargs["det_event"]


        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        if "tau" in kwargs:
            tau = kwargs["tau"]
        else:
            tau = 500
        if "tp" in kwargs:
            tp = kwargs["tp"]
        else:
            tp = 12

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=0, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # Exc pulse
                t=0, tp=tp, freq=0, flipangle=np.pi/2, 
                pcyc={"phases":[0, np.pi], "dets":[1, -1]}
            ))

        if hasattr(self, "pi_pulse"):
            pi_pulse = self.addPulse(self.pi_pulse.copy(
                t=tau, pcyc={"phases":[0], "dets": [1]}))
        else:
            pi_pulse = self.addPulse(RectPulse( # Pump 1 pulse
                t=tau, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*tau))
        else:
            self.addPulse(Detection(t=2*tau, tp=self.det_window.value))

# =============================================================================

class T1InversionRecoverySequence(Sequence):
    """
    Represents a T1 Inversion Recovery Sequence. 

    Sequence:
    [\pi - \tau - \pi/2 - \tau - echo]
    
    Parameters
    ----------
    B : int or float
        The B0 field, in Guass
    LO : int or float
        The LO frequency in GHz
    reptime : _type_
        The shot repetition time in us
    averages : int
        The number of scans.
    shots : int
        The number of shots per point
    start : float
        The minimum interpulse delay in ns, by default 300 ns
    step : float
        The step size of the interpulse delay in ns, by default 50 ns
    dim : int
        The number of points in the X axis

    Optional Parameters
    -------------------
    pi2_pulse : Pulse
        An autoEPR Pulse object describing the excitation pi/2 pulse. If
        not specified a RectPulse will be created instead. 
    pi_pulse : Pulse
        An autoEPR Pulse object describing the refocusing pi pulses. If
        not specified a RectPulse will be created instead. 
    """
# =============================================================================

class T2RelaxationSequence(HahnEchoSequence):
    """
    Represents a T2 relaxation sequence. A Hahn Echo where the interpulse delay increases
    
    Parameters
    ----------
    B : int or float
        The B0 field, in Guass
    freq : int or float
        The freq frequency in GHz
    reptime : _type_
        The shot repetition time in us
    averages : int
        The number of scans.
    shots : int
        The number of shots per point
    start : float
        The minimum interpulse delay in ns, by default 300 ns
    step : float
        The step size of the interpulse delay in ns, by default 50 ns
    dim : int
        The number of points in the X axis

    Optional Parameters
    -------------------
    pi2_pulse : Pulse
        An autoEPR Pulse object describing the excitation pi/2 pulse. If
        not specified a RectPulse will be created instead. 
    pi_pulse : Pulse
        An autoEPR Pulse object describing the refocusing pi pulses. If
        not specified a RectPulse will be created instead. 
    """

    def __init__(self, *, B, freq, reptime, averages, shots,start=500, step=40, dim=200, **kwargs) -> None:

        self.tau = Parameter(name="tau", value=start,step=step,dim=dim, unit="ns", description="The interpulse delay",virtual=True)
        super().__init__(B=B, freq=freq, reptime=reptime, averages=averages, shots=shots,tau=self.tau, **kwargs)

        self.name = "T2RelaxationSequence"
        self.evolution([self.tau])

    def simulate(self,ESEEM_depth=0.1, Tm=1e3):
        func = lambda x, a, tau, e: a*np.exp(-(x/tau)**e)
        xaxis = val_in_ns(self.tau)
        data = func(xaxis,1,Tm,1.6)
        data = add_phaseshift(data, 0.05)
        if ESEEM_depth != 0:
            data *= _gen_ESEEM(xaxis, 7.842, ESEEM_depth)
        return xaxis, data

# =============================================================================

class FieldSweepSequence(HahnEchoSequence):
    """
    Represents a Field Sweep (EDFS) sequence. 
    """
    def __init__(self, *, B, freq, Bwidth, reptime, averages, shots, **kwargs) -> None:
        """Build a Field Sweep (EDFS) sequence using either rectangular pulses or
        specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        Bwidth: int or float
            The width of the field sweep, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        
        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """
        super().__init__(
            B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.name = "FieldSweepSequence"


        self.B = Parameter(
            "B", value=B, start = -Bwidth/2, step=1, dim=Bwidth, unit="Gauss", description="Field sweep width"
        )
        
        self.evolution([self.B])
    
    def simulate(self):
        """
        Simulate the Field Sweep spectra as a Nitroxide radical.
        
        returns
        -------
        axis : np.ndarray
            The x-axis of the simulation
        data : np.ndarray
            The y-axis of the simulation

        """
        Vmodel = create_Nmodel(self.freq.value *1e3)
        axis = self.B.value + self.B.axis[0]['axis']
        sim_axis = axis * 0.1
        Boffset=0
        gy = 2.0061
        gz = 2.0021
        axy = 0.488
        az = 3.66
        GB = 0.45
        scale=1

        data = Vmodel(sim_axis,Boffset,gy,gz,axy,az,GB,scale)
        data = add_phaseshift(data, 0.05)
        return axis,data
        
# =============================================================================

class ReptimeScan(HahnEchoSequence):
    """
    Represents a reptime scan of a Hahn Echo Sequence. 
    """
    def __init__(self, *, B, freq, reptime, reptime_max, averages, shots, reptime_min=20, dim=100, **kwargs) -> None:
        """A Hahn echo sequence is perfomed with the shot repetition time increasing.1

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        freq : int or float
            The freq frequency in GHz
        reptime: float
            The default reptime, this is used for tuning pulses etc...
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        reptime_max : float
            The maximum shot repetition time in us
        reptime_min : float
            The minimum shot repetition time in us, by default 20 us
        dim : int
            The number of points in the reptime axis
        
        
        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """
        min_reptime = reptime_min
        dim = dim
        step  = (reptime_max-min_reptime)/dim
        step = np.around(step,decimals=-1)
        step = np.around(step,decimals=-1)
        reptime = Parameter(
            "reptime", reptime, start = min_reptime-reptime, step=step, dim=100, unit="us",
            description = "The shot repetition time")
        
        super().__init__(
            B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.name = "ReptimeScan"

        self.evolution([self.reptime])

    def simulate(self):
        """
        Simulates the reptime scan as an exponential recovery.

        .. math::
            V(t) = 1 - e^(-t/T1)

        Returns
        -------
        t : np.ndarray
            The x-axis of the simulation
        data : np.ndarray
            The y-axis of the simulation
        """
        def func(x,T1):
            return 1-np.exp(-x/T1)
        t = self.reptime.value + self.reptime.axis[0]['axis']
        T1 = 2000 #us

        data = func(t,T1)
        data = add_phaseshift(data, 0.05)
        return t, data

# =============================================================================

class CarrPurcellSequence(Sequence):
    """
    Represents a Carr-Purcell sequence. 
    """
    def __init__(self, *, B, freq, reptime, averages, shots,
             n,start=300,step=50, dim=100,**kwargs) -> None:
        """Build a Carr-Purcell dynamical decoupling sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        n : int
            The number refocusing pulses
        start : float
            The minimum interpulse delay in ns, by default 300 ns
        step : float
            The step size of the interpulse delay in ns, by default 50 ns
        dim : int
            The number of points in the X axis

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """

        name = "CarrPurcellSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.t = Parameter(name="tau", value=start,step=step,dim=dim, unit="ns",
            description="First interpulse delay", virtual=True)
        self.n = Parameter(name="n", value=n,
            description="The number of pi pulses", unit="None", virtual=True)
        self.dim = Parameter(name="dim", value=dim, unit="None",
            description="The number of points in the X axis", virtual=True)

        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]
        if "det_event" in kwargs:
            self.det_event = kwargs["det_event"]

        self._build_sequence()

    def _build_sequence(self):

        n = self.n.value
        # dt = 20
        # dim = np.floor((self.tau.value/(2*self.n.value) -deadtime)/dt)
        # # multipliers = [1]
        # # multipliers += [1+2*i for i in range(1,n)]
        # # multipliers += [2*n]

        # axis = np.arange(deadtime,tau/(2*n),10)

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=0, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # pi/2
                t=0, tp=16, freq=0, flipangle=np.pi/2,
                pcyc={"phases":[0, np.pi], "dets": [1, -1]}
            ))

        for i in range(n):
            if i==(n-1):
                phases = [0]
                dets = [1]
            elif i == (n-2):
                phases = [0, np.pi]
                dets = [1, 1]
            else:
                phases = [0, np.pi/2, np.pi, -np.pi/2]
                dets = [1,-1,1,-1]
            if hasattr(self, "pi_pulse"):
                self.addPulse(self.pi_pulse.copy(
                    t=self.t*(2*i + 1), pcyc={"phases":phases, "dets": dets}))
            else:
                self.addPulse(RectPulse(  # pi
                    t=self.t*(2*i + 1), tp=32, freq=0, flipangle=np.pi,
                    pcyc={"phases":phases, "dets": dets}
                ))
        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=self.t*(2*n)))
        else:
            self.addPulse(Detection(t=self.t*(2*n), tp=512))
        
        self.evolution([self.t])

    def simulate(self,T_CP=2.5e3, e=1.8):
        """
        Simulates the Carr-Purcell sequence as a single component stretched exponential decay.

        .. math::
            V(t) = e^(-(t/T_{CP})^e)
        
        Parameters
        ----------
        T_CP : float
            The Carr-Purcell relaxation time, by default 2.5e3 ns
        e : float
            The stretching exponent, by default 1.8
        
        Returns
        -------
        xaxis : np.ndarray
            The x-axis of the simulation
        data : np.ndarray
        """
        

        xaxis = val_in_ns(self.t)
        func = lambda x, a, tau, e: a*np.exp(-(x/tau)**e)
        data = func(xaxis,1,T_CP,e)
        data = add_phaseshift(data, 0.05)
        return xaxis, data
# =============================================================================

class ResonatorProfileSequence(Sequence):
    """
    Builds nutation based Resonator Profile sequence. 
    """

    def __init__(self, *, B, freq, reptime, averages, shots, fwidth=0.3,dtp=2, **kwargs) -> None:
        """Build a resonator profile nutation sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        Bwidth: int or float
            The width of the field sweep, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        fwidth: float
            The frequency width of the resonator profile in GHz, 
            by default 0.3GHz
        fstep: float
            The frequency step for the profile in GHz, by default 0.02GHz
            
        dtp: float
            The time step for the test pulse in ns, by default 2 ns

        Optional Parameters
        -------------------
        tau1: float
            The delay between the nutating pulse and the Hahn Echo, 
            by default 2000 ns
        tau2: float
            The interpulse delay in the Hahn Echo, 
            by default 500 ns
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """

        name = "ResonatorProfileSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.gyro = freq/B
        self.fwidth = Parameter('fwidth',fwidth,'GHz','Half the frequency sw')
        self.fstep = Parameter('fstep',kwargs.get('fstep',0.02),'GHz','Frequency step for the profile')
        self.dtp = Parameter('dtp',dtp,'ns','Time step for the pulse')

        self.kwargs = kwargs

        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]

        self._build_sequence()

    def _build_sequence(self):

        tau1 = self.kwargs.get("tau1",2000) #2000
        tau2 = self.kwargs.get("tau2",500) #500
        dim = np.floor(120/self.dtp.value).astype(int)
        tp = Parameter("tp", 0, step=self.dtp.value, dim=dim, unit="ns", description="Test Pulse length")
        fwidth= self.fwidth.value
        fstep = self.fstep.value
        dim = np.floor(fwidth*2/fstep)
        center_freq = self.freq.value
        self.freq = Parameter("freq", center_freq, start=-fwidth, step=fstep, dim=dim, unit="GHz", description="frequency")
        self.B = Parameter(
            "B",((center_freq)/self.gyro), start=-fwidth/self.gyro, step=fstep/self.gyro, dim=dim,
            unit="Guass",link=self.freq,description="B0 Field" )
        
        self.addPulse(RectPulse(  # Hard pulse
            t=0, tp=tp, freq=0, flipangle="Hard"
        ))

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=tau1, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # pi/2
            t=tau1, tp=16, freq=0, flipangle=np.pi/2, 
            pcyc={"phases":[0, np.pi], "dets": [1, -1]}
            ))
        
        if hasattr(self, "pi_pulse"):
            self.addPulse(self.pi_pulse.copy(
                t=tau1+tau2))
        else:
            self.addPulse(RectPulse(  # pi/2
            t=tau1+tau2, tp=32, freq=0, flipangle=np.pi
            ))

        self.addPulse(Detection(t=tau1+2*tau2, tp=64))


        self.pulses[0].scale.value = 1

        
        self.evolution([tp, self.freq])

    def simulate(self, Q=100, fc=None, nu1=75,damping=0.06):
        """
        Simulates a resonator profile sequence as a damped oscillation.

        Parameters
        ----------
        Q : int
            The quality factor of the resonator, by default 100
        fc : float
            The center frequency of the resonator, by default None. If None the freq frequency is used.
        nu1 : float
            The maximum amplitude of the resonator profile, by default 75 MHz. This is a linear frequency scale.
        damping : float
            The damping factor of the resonator, by default 0.06.
        
        """

        if fc is None:
            fc = self.freq.value

        xmin = self.freq.value - self.fwidth.value
        xmax = self.freq.value + self.fwidth.value


        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y

        mode = lambda x: lorenz_fcn(x, fc, fc/Q)
        axis = np.linspace(xmin,xmax)
        scale = nu1/mode(axis).max()
        self.mode = lambda x: lorenz_fcn(x, fc, fc/Q) * scale

        damped_oscilations = lambda x, f, c: np.cos(2*np.pi*f*x) * np.exp(-c*x)
        damped_oscilations_vec = np.vectorize(damped_oscilations)
        freq_axis = self.freq.value + self.freq.axis[0]['axis']
        freq_len = freq_axis.shape[0]
        tp_x = val_in_ns(self.pulses[0].tp)
        tp_len = tp_x.shape[0]
        nut_freqs = self.mode(freq_axis)

        damped_oscilations_vec
        data = damped_oscilations_vec(tp_x.reshape(tp_len,1),nut_freqs.reshape(1,freq_len)*1e-3,damping)
        return [tp_x, freq_axis], data

# =============================================================================

class TWTProfileSequence(Sequence):
    """
    Builds TWT based Resonator Profile sequence. 
    """
    
    def __init__(self,*,B,freq,reptime,averages=1,shots=100,**kwargs) -> None:

        name = "TWTProfileSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        self._build_sequence()

    def _build_sequence(self,):

        tau1=2000
        tau2=500

        self.addPulse(RectPulse(  # Hard pulse
            t=0, tp=4, freq=0, flipangle="Hard"
        ))

        self.addPulse(RectPulse(  # pi/2
            t=tau1, tp=16, freq=0, flipangle=np.pi/2
        ))
        self.addPulse(RectPulse(  # pi
            t=tau1+tau2, tp=32, freq=0, flipangle=np.pi
        ))

        self.addPulse(Detection(t=tau1+2*tau2, tp=512))


        self.pulses[0].scale.value = 1
        nut_axis = np.arange(0,66,2)
        self.addPulsesProg(
            [0],
            ["tp"],
            0,
            nut_axis)

        # Add amplitude sweep
        width= 0.3
        axis = np.arange(0,1.01,0.01)
        self.addPulsesProg(
                [0],
                ["scale"],
                1,
                axis)

# =============================================================================
