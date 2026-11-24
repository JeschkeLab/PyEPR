from pyepr.classes import Parameter
from pyepr.pulses import Pulse, Detection, Delay, RectPulse
from pyepr.utils import val_in_ns, val_in_us, add_phaseshift, _gen_ESEEM
from pyepr.fieldsweep_analysis import create_Nmodel
from pyepr.plotting import draw_sequence
import pyepr.pulses as ad_pulses
import numpy as np
import matplotlib.pyplot as plt
import json
from itertools import product
import copy
from pyepr.utils import build_table, autoEPRDecoder
import json
from pyepr._version import __version__
import uuid
import base64
import numbers
import warnings
from typing import Callable

__all__ = [
    'Sequence', 'HahnEchoSequence', 'T1InversionRecoverySequence',
    'HahnEchoRelaxationSequence', 'T2RelaxationSequence', 'FieldSweepSequence',
    'ReptimeScanSequence', 'ReptimeScan', 'CarrPurcellSequence', 'ResonatorProfileSequence',
    'TWTProfileSequence', 'KBBSequence', 'KBBPulseProfileSequence',
]

class Sequence:
    """
    Represents an experimental pulse sequence.
    """

    time: Parameter
    """The estimated run time of the sequence (HH:MM:SS), updated whenever the sequence changes."""
    pcyc_vars: list[int]
    """The indices of the pulses that are phase cycled."""
    pcyc_cycles: np.ndarray
    """The phase of each phase-cycled pulse for every step of the phase cycle."""
    pcyc_dets: np.ndarray
    """The detection sign for every step of the phase cycle."""

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

        Examples
        --------
        A custom sequence can be built by adding pulses to an empty sequence.

        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np
            epr.set_waveform_precision(1)

            seq = epr.Sequence(
                name="CustomHahnEcho", B=12000, freq=34, reptime=1e3,
                averages=1, shots=10)
            seq.addPulse(epr.RectPulse(t=0, tp=16, freq=0, flipangle=np.pi/2))
            seq.addPulse(epr.RectPulse(t=500, tp=32, freq=0, flipangle=np.pi))
            seq.addPulse(epr.Detection(t=1000, tp=512))
            seq.draw()
        """

        self.pulses: list[Pulse] = []
        """The list of pulses in this sequence."""

        self.uuid = uuid.uuid1()
        """Unique identifier of this sequence. Reset when the sequence is
        copied, a pulse is added or the evolution is changed."""

        self._axes_uuid = []
        self._reduce_uuid = []

        if isinstance(B, Parameter):
            self.B: Parameter = B.copy()
            """The static B0 field for the experiment in Gauss."""
        else:
            self.B: Parameter = Parameter(
                "B", B, "Gauss",
                "The static B0 field for the experiment")
            """The static B0 field for the experiment in Gauss."""

        if freq is not None:
            self.freq: Parameter = Parameter(
                "freq", freq, "GHz",
                "The central frequency of the sequence")
            """The central frequency of the sequence in GHz."""
        elif LO is not None:
            self.freq: Parameter = Parameter(
                "freq", LO, "GHz",
                "The central frequency of the sequence")
            """The central frequency of the sequence in GHz."""
            warnings.warn("LO is deprecated, please use freq instead", DeprecationWarning, stacklevel=2)
        else:
            raise ValueError("Either freq or LO must be specified")

        if isinstance(reptime, Parameter):
            self.reptime: Parameter = reptime.copy()
            """The shot repetition time in us."""
        else:
            self.reptime: Parameter = Parameter(
                "reptime", reptime, "us",
                "The shot repetition time")
            """The shot repetition time in us."""

        self.averages: Parameter = Parameter("averages", int(averages), "None",
            "The number of averages to perform.")
        """The number of averages (scans) to perform."""

        self.shots: Parameter = Parameter(
            "shots", int(shots), "None",
            "The number of shots per scan.")
        """The number of shots per point."""

        self.name: str = name
        """The name of this pulse sequence."""
        self.progTable: dict = {"EventID": [], "Variable": [], "axis": [],
                     "axID": [], "uuid": [], "reduce": []}
        """The program table for this sequence, which is used to generate the pulse program. Has columns for EventID, Variable, axis, axID, uuid, and reduce."""
        self.evo_params: list[Parameter] = []
        """The list of evolution parameters for this sequence, which are used to define the evolution axes of the sequence."""
        
        pass


    @property
    def num_pulses(self) -> int:
        """The number of pulses in this sequence."""
        return len(self.pulses)
    
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

    def draw(self, **kwargs):
        """
        Draws the sequence as an annotated matplotlib figure, showing pictorial representations of the pulses and delays in the sequence.
        Calls the `draw_sequence` function from the `pyepr.plotting` module.
        """
        return draw_sequence(self, **kwargs)

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
        self.uuid = uuid.uuid1() # UUID reset when the pulses change
        self._buildPhaseCycle()
        self._estimate_time()

    def _estimate_time(self) -> float:
        """
        Calculates the estimated experiment time in seconds.
        """
        self._buildPhaseCycle()
        acqs = self.averages.value * self.shots.value * self.pcyc_dets.shape[0]
        
        reptime = np.mean(val_in_us(self.reptime))
        if hasattr(self,'evo_params'):
            acqs *= np.prod([np.prod(param.dim) for param in self.evo_params])        
        time = acqs * reptime * 1e-6

        self.time = Parameter(name="time", value=f"{(time // 3600):.0f}:{(time % 3600) // 60:.0f}:{(time % 60):.0f}", unit="HH:MM:SS",
        description="Estimated sequence run time")
        return time

    def _buildPhaseCycle(self):
        # Group phase-cycled pulses by link target, so linked pulses step
        # through their phase cycle in phase (zipped) instead of being
        # combined via a full outer product with the rest of the sequence.
        groups = {}
        for pulse in self.pulses:
            if pulse.pcyc is None:
                continue
            link = pulse.pcyc.get("link")
            key = pulse.uuid if link is None else link
            groups.setdefault(key, []).append(pulse)
        groups = list(groups.values())

        phase_axes = [list(zip(*(p.pcyc["Phases"] for p in g))) for g in groups]
        det_axes = [list(zip(*(p.pcyc["DetSigns"] for p in g))) for g in groups]

        self.pcyc_vars = [self.pulses.index(p) for g in groups for p in g]
        self.pcyc_cycles = np.array([
            [phase for phases in row for phase in phases]
            for row in product(*phase_axes)])
        self.pcyc_dets = np.array([
            np.prod([det for dets in row for det in dets])
            for row in product(*det_axes)])

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
        self.uuid = uuid.uuid1() # UUID reset when the evolution changes
        self._axes_uuid = [param.uuid for param in params]
        self._reduce_uuid = [param.uuid for param in reduce]

        self._buildProgTable()


    def _buildProgTable(self):

        progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": [], "uuid": [], "reduce": []}
        
        for n, pulse in enumerate(self.pulses):
            table = pulse.build_table()
            for i in range(len(table["uuid"])):
                if table["uuid"][i] in self._axes_uuid:
                    progTable["axID"].append(self._axes_uuid.index(table["uuid"][i]))
                    progTable["uuid"].append(table["uuid"][i]) 
                    progTable["EventID"].append(n)
                    progTable["Variable"].append(table["Variable"][i])
                    progTable["axis"].append(table["axis"][i])
                    if table["uuid"][i] in self._reduce_uuid:
                        progTable["reduce"].append(True)
                    else:
                        progTable["reduce"].append(False)

        for var_name in vars(self):
            var = getattr(self, var_name)
            if type(var) is Parameter:
                if not var.is_static() and not var.virtual:
                    for i in range(len(var.axis)):
                        if var.axis[i]["uuid"] in self._axes_uuid:
                            progTable["axID"].append(self._axes_uuid.index(var.axis[i]["uuid"]))
                            progTable["EventID"].append(None)
                            progTable["Variable"].append(var_name) 
                            progTable["axis" ].append(var.axis[i]["axis"])
                            progTable["uuid"].append(var.axis[i]["uuid"]) 
                            if var.axis[i]["uuid"] in self._reduce_uuid:
                                progTable["reduce"].append(True)
                            else:
                                progTable["reduce"].append(False)
        self.progTable = progTable
        self._estimate_time()
        return self.progTable
    
    @property
    def seqtable_steps(self) -> int:
        """The number of steps in the sequence table, which is the product of the number of phase cycle detections, the number of phase cycle variables, and the product of the dimensions of the evolution parameters."""
        if len(self.evo_params) > 0:
            return self.pcyc_dets.shape[0] * (len(self.pcyc_vars)+1) * np.prod([np.prod(param.dim) for param in self.evo_params])
    
    @property
    def shape(self) -> list[int]:
        """
        Gives the shape of the sequence,excluding the number of shots.

        Returns
        -------
        list
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

        # Flip the axes dim so it is Z,Y,X...
        axes_dim = axes_dim[::-1]

        shape = [nAvgs]+  axes_dim +[nPcyc]

        # Force list output to all be int
        shape = np.array(shape).astype(int).tolist()

        return shape



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
        """Shifts the detection frequency of all pulses to zero, by subtracting the detection frequency from all pulses and the sequence frequency."""
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
                elif param.value is None:
                    seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
                        param.name, "None", unit, param.description)
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

    def copy(self, reset_uuid=True):
        """Deep copy the sequence object, including all pulses and parameters.

        Parameters
        ----------
        reset_uuid : bool, optional
            If True, the UUID of the copied sequence will be reset to a new
            value. Default is True.
        """
        new_sequence = copy.deepcopy(self)
        if reset_uuid:
            new_sequence.uuid = uuid.uuid1()
        return new_sequence
    
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
        r"""Build a Hahn-Echo sequence using either rectangular pulses or
        specified pulses. By default no progression is added to this sequence.

        Sequence:
        [\pi/2 - \tau - \pi - \tau - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
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
        tau : float or Parameter
            The interpulse delay in ns, by default 500 ns
        tp : float
            The length of the default rectangular pulses in ns, by default 12 ns
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

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np

            epr.set_waveform_precision(1)
            seq = epr.HahnEchoSequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, tau=500)
            seq.draw()


        """

        name = kwargs.pop('name',"HahnEchoSequence")

        if "pi_pulse" in kwargs:
            self.pi_pulse: Pulse = kwargs["pi_pulse"]
            """The refocusing pi pulse."""
        if "pi2_pulse" in kwargs:
            self.pi2_pulse: Pulse = kwargs["pi2_pulse"]
            """The excitation pi/2 pulse."""
        if "det_event" in kwargs:
            self.det_event: Detection = kwargs["det_event"]
            """The detection event."""


        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        tau = kwargs.get("tau", 500)
        if isinstance(tau, Parameter):
            self.tau: Parameter = tau
            """The interpulse delay in ns."""
        else:
            self.tau: Parameter = Parameter(
                name="tau", value=tau, unit="ns",
                description="The interpulse delay", virtual=True)
            """The interpulse delay in ns."""
        tp = kwargs.get("tp", 12)

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
                t=self.tau, pcyc={"phases":[0], "dets": [1]}))
        else:
            pi_pulse = self.addPulse(RectPulse( # Pump 1 pulse
                t=self.tau, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*self.tau))
        else:
            self.addPulse(Detection(t=2*self.tau, tp=512))

    def simulate(self, Tm=1e3):
        """
        Simulates the Hahn echo intensity as a stretched exponential decay
        of the interpulse delay, with a stretching exponent of 1.6.

        .. math::
            V(\\tau) = e^{-(\\tau/T_m)^{1.6}}

        Parameters
        ----------
        Tm : float
            The phase memory time in ns, by default 1e3 ns

        Returns
        -------
        xaxis : np.ndarray
            The interpulse delay in ns (a single point)
        data : np.ndarray
            The simulated complex echo intensity (a single point)

        Examples
        --------
        .. code-block:: python

            import pyepr as epr
            seq = epr.HahnEchoSequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, tau=500)
            tau, data = seq.simulate(Tm=1e3)
            print(tau, data)
        """
        xaxis = np.atleast_1d(val_in_ns(self.tau))
        data = np.exp(-(xaxis/Tm)**1.6)
        data = add_phaseshift(data, 0.05)
        return xaxis, data

# =============================================================================

class T1InversionRecoverySequence(Sequence):
    """
    Represents a T1 Inversion Recovery Sequence.
    """

    def __init__(self, *, B, freq, reptime, averages, shots,start=500, step=40, dim=200, **kwargs) -> None:
        r"""Build a T1 inversion recovery sequence using either rectangular
        pulses or specified pulses.

        Sequence:
        [\pi - t - \pi/2 - \tau - \pi - \tau - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        start : float
            The minimum inversion recovery delay in ns, by default 500 ns
        step : float
            The step size of the inversion recovery delay in ns, by default 40 ns
        dim : int
            The number of points in the X axis, by default 200

        Optional Parameters
        -------------------
        tau : float
            The Hahn echo interpulse delay in ns, by default 500 ns
        tp : float
            The length of the default rectangular pulses in ns, by default 12 ns
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead.
        pi_pulse : Pulse
            An autoEPR Pulse object describing the inversion and refocusing
            pi pulses. If not specified a RectPulse will be created instead.
        det_event : Pulse
            An autoEPR Pulse object describing the detection event. If
            not specified a standard detection event will be created instead.

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)
            seq = epr.T1InversionRecoverySequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=10000, step=2500, dim=500)
            seq.draw()

        """
        name = "T1InversionRecoverySequence"
        tp = kwargs.get("tp", 12)

        self.tau: Parameter = Parameter(name="tau", value=kwargs.get("tau", 500), unit="ns", description="The hahn echo interpulse delay",virtual=True)
        """The Hahn echo interpulse delay in ns."""
        self.t: Parameter = Parameter(name="t", value=start, step=step, dim=dim, unit="ns", description="The inversion recovery delay",virtual=True)
        """The inversion recovery delay in ns."""
        super().__init__(name=name,B=B, freq=freq, reptime=reptime, averages=averages, shots=shots, **kwargs)

        if "pi_pulse" in kwargs:
            self.pi_pulse: Pulse = kwargs["pi_pulse"]
            """The inversion and refocusing pi pulse."""
        if "pi2_pulse" in kwargs:
            self.pi2_pulse: Pulse = kwargs["pi2_pulse"]
            """The excitation pi/2 pulse."""
        if "det_event" in kwargs:
            self.det_event: Detection = kwargs["det_event"]
            """The detection event."""

        if hasattr(self, "pi_pulse"):
            pi_pulse = self.addPulse(self.pi_pulse.copy(
                t=0, pcyc={"phases":[0], "dets": [1]}))
        else:
            pi_pulse = self.addPulse(RectPulse( # Pump 1 pulse
                t=0, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=self.t, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # Exc pulse
                t=self.t, tp=tp, freq=0, flipangle=np.pi/2, 
                pcyc={"phases":[0, np.pi], "dets":[1, -1]}
            ))

        if hasattr(self, "pi_pulse"):
            pi_pulse = self.addPulse(self.pi_pulse.copy(
                t=self.t+self.tau, pcyc={"phases":[0], "dets": [1]}))
        else:
            pi_pulse = self.addPulse(RectPulse( # Pump 1 pulse
                t=self.t+self.tau, tp=tp, freq=0, flipangle=np.pi
            ))
        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=self.t+2*self.tau))
        else:
            self.addPulse(Detection(t=self.t+2*self.tau, tp=self.pulses[1].tp.value))

        self.evolution([self.t])

    def simulate(self, T1=1e4):
        """
        Simulates the T1 recovery as a simple exponential recovery.

        .. math::
            V(t) = 1 - 2e^{-t/T_1}

        Parameters
        ----------
        T1 : float
            The longitudinal relaxation time in ns, by default 1e4 ns

        Returns
        -------
        xaxis : np.ndarray
            The inversion recovery delay in ns
        data : np.ndarray
            The simulated complex echo intensity

        Examples
        --------
        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.T1InversionRecoverySequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=500, step=500, dim=100)
            x, data = seq.simulate(T1=1e4)
            plt.plot(x, data.real)
            plt.xlabel("Inversion Recovery Delay (ns)")
            plt.ylabel("Signal Intensity (a.u.)")
        """
        func = lambda x, a, T1: a*(1 - 2*np.exp(-x/T1))
        xaxis = val_in_ns(self.t)
        data = func(xaxis,1,T1)
        data = add_phaseshift(data, 0.05)
        return xaxis, data
# =============================================================================


class HahnEchoRelaxationSequence(HahnEchoSequence):
    """
    Represents a Hahn Echo relaxation sequence for measureing Tm. A Hahn Echo where the interpulse delay increases
    """

    def __init__(self, *, B, freq, reptime, averages, shots,start=500, step=40, dim=200, **kwargs) -> None:
        r"""Build a Hahn echo relaxation sequence, where the interpulse delay
        is incremented, using either rectangular pulses or specified pulses.

        Sequence:
        [\pi/2 - \tau - \pi - \tau - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        start : float
            The minimum interpulse delay in ns, by default 500 ns
        step : float
            The step size of the interpulse delay in ns, by default 40 ns
        dim : int
            The number of points in the X axis, by default 200

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead.
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead.

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)
            seq = epr.HahnEchoRelaxationSequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=300, step=40, dim=200)
            seq.draw()

        """

        self.tau: Parameter = Parameter(name="tau", value=start,step=step,dim=dim, unit="ns", description="The interpulse delay",virtual=True)
        """The incremented interpulse delay in ns."""
        super().__init__(B=B, freq=freq, reptime=reptime, averages=averages, shots=shots,tau=self.tau, **kwargs)

        self.name = "HahnEchoRelaxationSequence"
        self.evolution([self.tau])

    def simulate(self,ESEEM_depth=0.1, Tm=1e3):
        """
        Simulates the Hahn echo decay as a stretched exponential, with a
        stretching exponent of 1.6, modulated by proton ESEEM.

        .. math::
            V(\\tau) = e^{-(\\tau/T_m)^{1.6}}

        Parameters
        ----------
        ESEEM_depth : float
            The depth of the proton ESEEM modulation, by default 0.1. Set to 0
            to disable the ESEEM modulation.
        Tm : float
            The phase memory time in ns, by default 1e3 ns

        Returns
        -------
        xaxis : np.ndarray
            The interpulse delay in ns
        data : np.ndarray
            The simulated complex echo intensity

        Examples
        --------
        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.HahnEchoRelaxationSequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=300, step=40, dim=200)
            x, data = seq.simulate(ESEEM_depth=0.1, Tm=1e3)
            plt.plot(x, data.real)
            plt.xlabel("Interpulse Delay (ns)")
            plt.ylabel("Signal Intensity (a.u.)")
        """
        func = lambda x, a, tau, e: a*np.exp(-(x/tau)**e)
        xaxis = val_in_ns(self.tau)
        data = func(xaxis,1,Tm,1.6)
        data = add_phaseshift(data, 0.05)
        if ESEEM_depth != 0:
            data *= _gen_ESEEM(xaxis, 7.842, ESEEM_depth)
        return xaxis, data

class T2RelaxationSequence(HahnEchoRelaxationSequence):
    """
    Represents a T2 relaxation sequence. This is an alias for HahnEchoRelaxationSequence.
    A Hahn Echo where the interpulse delay increases.
    """

    def __init__(self, *, B, freq, reptime, averages, shots, start=500, step=40, dim=200, **kwargs) -> None:
        r"""Build a T2 relaxation sequence. This is an alias for
        :class:`HahnEchoRelaxationSequence`.

        Sequence:
        [\pi/2 - \tau - \pi - \tau - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        start : float
            The minimum interpulse delay in ns, by default 500 ns
        step : float
            The step size of the interpulse delay in ns, by default 40 ns
        dim : int
            The number of points in the X axis, by default 200

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead.
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead.

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)
            seq = epr.T2RelaxationSequence(
                B=12000, freq=34, reptime=1e3, averages=1, shots=10, start=300, step=40, dim=200)
            seq.draw()

        """
        super().__init__(
            B=B, freq=freq, reptime=reptime, averages=averages, shots=shots,
            start=start, step=step, dim=dim, **kwargs)

# =============================================================================

class FieldSweepSequence(HahnEchoSequence):
    """
    Represents a Field Sweep (EDFS) sequence. 

    

    """
    def __init__(self, *, B, freq, Bwidth, reptime, averages, shots, **kwargs) -> None:
        r"""Build a Field Sweep (EDFS) sequence using either rectangular pulses or
        specified pulses.

        Sequence:
        [\pi/2 - \tau - \pi - \tau - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
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
        step : float
            The field step in Gauss, by default 1 G
        tau : float
            The Hahn echo interpulse delay in ns, by default 500 ns
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead.
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead.

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.FieldSweepSequence(
            B=12000, Bwidth=100, freq=34, reptime=1e3, averages=1, shots=10)
            seq.draw()

        """
        super().__init__(
            B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.name = "FieldSweepSequence"

        dim = Bwidth // kwargs.get('step',1)
        self.B: Parameter = Parameter(
            "B", value=B, start = -Bwidth/2, step=kwargs.get('step',1), dim=dim, unit="Gauss", description="Field sweep width"
        )
        """The swept B0 field in Gauss, centred on the given field."""

        self.evolution([self.B])

    def simulate(self):
        """
        Simulate the Field Sweep spectra as a Nitroxide radical.

        Returns
        -------
        axis : np.ndarray
            The magnetic field axis in Gauss
        data : np.ndarray
            The simulated complex echo intensity

        Examples
        --------
        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            epr.set_waveform_precision(1)

            N_gyro = 0.002803632236095
            freq = 34
            seq = epr.FieldSweepSequence(
                B=freq/N_gyro, Bwidth=200, freq=freq, reptime=1e3, averages=1, shots=10)
            B, data = seq.simulate()
            plt.plot(B, data.real)
            plt.xlabel("Magnetic Field (G)")
            plt.ylabel("Signal Intensity (a.u.)")
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


class ReptimeScanSequence(HahnEchoSequence):
    """
    Represents a reptime scan of a Hahn Echo Sequence. 
    """
    def __init__(self, *, B, freq, reptime, reptime_max, averages, shots, start=20, dim=100, **kwargs) -> None:
        r"""A Hahn echo sequence is perfomed with the shot repetition time increasing.1

        Sequence:
        [\pi/2 - \tau - \pi - \tau - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
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
        start : float
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

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)
            seq = epr.ReptimeScanSequence(
            B=12000, freq=34, reptime=1000, averages=1,
            shots=10, reptime_max=10000, start=20, dim=100)  

            seq.draw()
        """
        min_reptime = start
        step  = (reptime_max-min_reptime)/dim
        step = np.around(step,decimals=-1)
        step = np.around(step,decimals=-1)
        reptime = Parameter(
            "reptime", reptime, start = min_reptime-reptime, step=step, dim=dim, unit="us",
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
            The shot repetition time axis in us
        data : np.ndarray
            The simulated complex echo intensity

        Examples
        --------
        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.ReptimeScanSequence(
                B=12000, freq=34, reptime=1000, averages=1,
                shots=10, reptime_max=10000, start=20, dim=100)
            t, data = seq.simulate()
            plt.plot(t, data.real)
            plt.xlabel(r"Shot Repetition Time ($\\mu$s)")
            plt.ylabel("Signal Intensity (a.u.)")
        """
        def func(x,T1):
            return 1-np.exp(-x/T1)
        t = self.reptime.value + self.reptime.axis[0]['axis']
        T1 = 2000 #us

        data = func(t,T1)
        data = add_phaseshift(data, 0.05)
        return t, data

class ReptimeScan(ReptimeScanSequence):
    """
    Deprecated alias for :class:`ReptimeScanSequence`.
    """

    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "ReptimeScan"
        print("Warning: ReptimeScan is deprecated, use ReptimeScanSequence instead.")
# =============================================================================

class CarrPurcellSequence(Sequence):
    """
    Represents a Carr-Purcell sequence. 
    """
    def __init__(self, *, B, freq, reptime, averages, shots,
             n,start=300,step=50, dim=100,**kwargs) -> None:
        r"""Build a Carr-Purcell dynamical decoupling sequence using either 
        rectangular pulses or specified pulses.

        Sequence:
        [\pi/2 [- \tau - \pi - \tau]^n - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
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

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            epr.set_waveform_precision(1)
            seq = epr.CarrPurcellSequence(B=12000, freq=34, reptime=1000, averages=1, shots=10, n=4)  

            seq.draw()
        """

        name = "CarrPurcellSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.t: Parameter = Parameter(name="tau", value=start,step=step,dim=dim, unit="ns",
            description="First interpulse delay", virtual=True)
        """The first interpulse delay in ns, which is incremented."""
        self.n: Parameter = Parameter(name="n", value=n,
            description="The number of pi pulses", unit="None", virtual=True)
        """The number of refocusing pi pulses."""
        self.dim: Parameter = Parameter(name="dim", value=dim, unit="None",
            description="The number of points in the X axis", virtual=True)
        """The number of points in the X axis."""

        if "pi_pulse" in kwargs:
            self.pi_pulse: Pulse = kwargs["pi_pulse"]
            """The refocusing pi pulse."""
        if "pi2_pulse" in kwargs:
            self.pi2_pulse: Pulse = kwargs["pi2_pulse"]
            """The excitation pi/2 pulse."""
        if "det_event" in kwargs:
            self.det_event: Detection = kwargs["det_event"]
            """The detection event."""

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
            V(t) = e^{(-(t/T_{CP})^e)}

        Examples
        --------

        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.CarrPurcellSequence(B=12000, freq=34, reptime=1000, averages=1, shots=10, n=4)  
            x,data = seq.simulate(T_CP=2.5e3, e=1.8)
            plt.plot(x,data)
            plt.xlabel("Interpulse Delay (ns)")
            plt.ylabel("Signal Intensity (a.u.)")    

        
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
        r"""Build a resonator profile nutation sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
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
        step:
            The frequency step parameter in GHz, by default 0.02GHz
        dim: int
            The number of points in the tp axis, by default set to 120/dtp, which is 60 points for dtp=2ns.

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

        Examples
        --------

        .. plot::
            :include-source:
            
            import pyepr as epr
            import numpy as np
            epr.set_waveform_precision(1)

            N_gyro = 0.002803
            freq = 34
            pi2_pulse = epr.RectPulse(tp=32,freq=0,flipangle=np.pi/2)
            pi_pulse = epr.RectPulse(tp=32,freq=0,flipangle=np.pi)

            seq = epr.ResonatorProfileSequence(
                B=freq/N_gyro,
                freq=freq,
                reptime=3e3,
                averages=1,
                shots=10,
                pi2_pulse=pi2_pulse,
                pi_pulse=pi_pulse,
                fwidth=0.3,
                fstep=0.02,
                dtp=1,
                dim=40)

            seq.draw()

        """

        name = "ResonatorProfileSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.gyro: float = freq/B
        """Gyromagenetic ratio in GHz/Gauss"""
        self.fwidth: Parameter = Parameter('fwidth',fwidth,'GHz','Half the frequency sw')
        """Half the frequency sweep width for the profile"""
        self.fstep: Parameter = Parameter('fstep',kwargs.get('fstep',0.02),'GHz','Frequency step for the profile')
        """Frequency step for the profile"""
        self.dtp: Parameter = Parameter('dtp',float(dtp),'ns','Time step for the pulse')
        """Time step for the test pulse"""
        self.kwargs: dict = kwargs
        """Keyword arguments"""

        if "pi_pulse" in kwargs:
            self.pi_pulse: Pulse = kwargs["pi_pulse"]
            """The refocusing pi pulse of the Hahn echo detection."""
        if "pi2_pulse" in kwargs:
            self.pi2_pulse: Pulse = kwargs["pi2_pulse"]
            """The excitation pi/2 pulse of the Hahn echo detection."""

        self._build_sequence()

    def _build_sequence(self):

        tau1 = self.kwargs.get("tau1",2000) #2000
        tau2 = self.kwargs.get("tau2",500) #500
        fstep = self.kwargs.get("step",0.02)
        dim = self.kwargs.get("dim",np.floor(120/self.dtp.value).astype(int))
        tp = Parameter("tp", 0, step=self.dtp.value, dim=dim, unit="ns", description="Test Pulse length")
        fwidth= self.fwidth.value
        fstep = self.fstep.value
        dim = np.floor(fwidth*2/fstep)
        center_freq = self.freq.value
        self.freq: Parameter = Parameter("freq", center_freq, start=-fwidth, step=fstep, dim=dim, unit="GHz", description="frequency")
        """The swept central frequency of the sequence in GHz."""
        self.B: Parameter = Parameter(
            "B",((center_freq)/self.gyro), start=-fwidth/self.gyro, step=fstep/self.gyro, dim=dim,
            unit="Gauss",link=self.freq,description="B0 Field" )
        """The B0 field in Gauss, swept together with the frequency."""

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

        self.addPulse(Detection(t=tau1+2*tau2, tp=512))


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

        Returns
        -------
        axes : list[np.ndarray]
            The test pulse length axis in ns and the frequency axis in GHz
        data : np.ndarray
            The simulated nutation data, with shape (len(tp), len(freq))

        Examples
        --------
        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            import numpy as np
            epr.set_waveform_precision(1)

            N_gyro = 0.002803
            freq = 34
            seq = epr.ResonatorProfileSequence(
                B=freq/N_gyro, freq=freq, reptime=3e3, averages=1, shots=10,
                fwidth=0.3, fstep=0.02, dtp=1, dim=100)
            [tp, f], data = seq.simulate(Q=100, nu1=75)
            plt.pcolormesh(f, tp, data.real, shading='auto')
            plt.xlabel("Frequency (GHz)")
            plt.ylabel("Test Pulse Length (ns)")
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
        self.mode: Callable[[np.ndarray], np.ndarray] = lambda x: lorenz_fcn(x, fc, fc/Q) * scale
        """The simulated resonator mode, giving the nutation frequency in MHz as a function of frequency in GHz."""

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
    
    def __init__(self,*,B,freq,reptime,averages=1,shots=100,dtp=2,**kwargs) -> None:
        r"""Build a TWT profile sequence, a nutation experiment where both the
        length and the amplitude scale of the nutating pulse are varied. This
        is used to characterise the (non-)linearity of the TWT amplifier.

        Sequence:
        [p_{nut}(t_p, scale) - \tau_1 - \pi/2 - \tau_2 - \pi - \tau_2 - echo]

        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int, optional
            The number of scans, by default 1
        shots : int, optional
            The number of shots per point, by default 100
        dtp : float, optional
            The time step for the test pulse in ns, by default 2 ns

        Optional Parameters
        -------------------
        tau1 : float
            The delay between the nutating pulse and the Hahn Echo,
            by default 2000 ns
        tau2 : float
            The interpulse delay in the Hahn Echo, by default 500 ns
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead.
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead.

        Examples
        --------
        .. plot::
            :include-source:

            import pyepr as epr
            import numpy as np
            epr.set_waveform_precision(1)

            seq = epr.TWTProfileSequence(
                B=12000, freq=34, reptime=3e3, averages=1, shots=10, dtp=2)
            seq.draw()

        """

        name = "TWTProfileSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        self.kwargs: dict = kwargs
        """Keyword arguments"""
        self.dtp: Parameter = Parameter('dtp',dtp,'ns','Time step for the pulse')
        """Time step for the test pulse"""

        if "pi_pulse" in kwargs:
            self.pi_pulse: Pulse = kwargs["pi_pulse"]
            """The refocusing pi pulse of the Hahn echo detection."""
        if "pi2_pulse" in kwargs:
            self.pi2_pulse: Pulse = kwargs["pi2_pulse"]
            """The excitation pi/2 pulse of the Hahn echo detection."""

        self._build_sequence()

    def _build_sequence(self,):

        tau1 = self.kwargs.get("tau1",2000) #2000
        tau2 = self.kwargs.get("tau2",500) #500

        dim = np.floor(120/self.dtp.value).astype(int)
        tp = Parameter("tp", 0, step=self.dtp.value, dim=dim, unit="ns", description="Test Pulse length")
        scale = Parameter("scale", 0, step=0.01, dim=100, unit="None", description="Amplitude scale factor")

        self.addPulse(RectPulse(  # Hard pulse
            t=0, tp=tp, freq=0, flipangle="Hard",scale=scale
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

        self.addPulse(Detection(t=tau1+2*tau2, tp=512))

        self.evolution([tp, scale])

    def simulate(self, nu1=75, damping=0.06, sat=None):
        """
        Simulates the TWT profile as damped nutations, where the nutation
        frequency is proportional to the amplitude scale of the test pulse.

        An optional saturation of the amplifier can be included, in which case
        the nutation frequency follows a hyperbolic tangent of the scale.

        .. math::
            \\nu_1(s) = \\nu_{1,max} \\frac{\\tanh(s/s_{sat})}{\\tanh(1/s_{sat})}

        Parameters
        ----------
        nu1 : float
            The nutation frequency at full scale (scale=1) in MHz, by default 75 MHz
        damping : float
            The damping factor of the nutations, by default 0.06
        sat : float, optional
            The saturation scale of the amplifier. If None the amplifier is
            assumed to be linear, by default None

        Returns
        -------
        axes : list[np.ndarray]
            The test pulse length axis in ns and the amplitude scale axis
        data : np.ndarray
            The simulated nutation data, with shape (len(tp), len(scale))

        Examples
        --------
        .. plot::
            :include-source:

            import matplotlib.pyplot as plt
            import pyepr as epr
            epr.set_waveform_precision(1)

            seq = epr.TWTProfileSequence(
                B=12000, freq=34, reptime=3e3, averages=1, shots=10, dtp=2)
            [tp, scale], data = seq.simulate(nu1=75, sat=0.5)
            plt.pcolormesh(scale, tp, data.real, shading='auto')
            plt.xlabel("Amplitude Scale")
            plt.ylabel("Test Pulse Length (ns)")
        """
        tp_x = val_in_ns(self.pulses[0].tp)
        scale_x = self.pulses[0].scale.value + self.pulses[0].scale.axis[0]['axis']

        if sat is None:
            nut_freqs = nu1 * scale_x
        else:
            nut_freqs = nu1 * np.tanh(scale_x/sat) / np.tanh(1/sat)

        data = np.cos(2*np.pi*tp_x.reshape(-1,1)*nut_freqs.reshape(1,-1)*1e-3) \
            * np.exp(-damping*tp_x.reshape(-1,1))
        data = add_phaseshift(data, 0.05)
        return [tp_x, scale_x], data

# =============================================================================

class KBBSequence(Sequence):
    """
    Builds a Kutnz-Boehlen-Bodenhausen echo sequence.

    """

    def __init__(self,*,B,freq,reptime,averages=1,shots=100,**kwargs) -> None:
        """
        Parameters
        ----------
        B : int or float
            The B0 field, in Gauss
        freq : int or float
            The freq frequency in GHz
        reptime : float
            The shot repetition time in us
        averages : int, optional
            The number of scans, by default 1
        shots : int, optional
            The number of shots per point, by default 100
        """
        name = kwargs.pop("name", "KBBSequence")
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        self.kwargs: dict = kwargs
        """Keyword arguments"""

        self._build_sequence()

    def _build_sequence(self,):
        pass

# =============================================================================

class KBBPulseProfileSequence(KBBSequence):
    """
    Builds a KBB based Resonator Profile sequence. 
    """
    
    def __init__(self,*,B,freq,reptime,averages=1,shots=100,dtp=2,**kwargs) -> None:

        name = "KBB_PulseProfileSequence"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self._build_sequence()

    def _build_sequence(self,):
        pass