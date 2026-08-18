from pyepr import Sequence, Parameter, Detection, ChirpPulse, RectPulse

class LoSweepTuneSequence(Sequence):

    def __init__(self,*,freq, **kwargs):
        """
        LO Sweep Tune Sequence

        This is formed of a single monochromatic rectangular pulse which is detected, at different LO frequencies.
        The resonator dip is charecterised by the absorption of the pulse, and less reflected power.
        
        Parameters
        ----------
        freq : float
            The center frequency in GHz.
        
        B: float
            The magnetic field in gAUSS, default is 0.
        reptime : float
            The repetition time in us, default is 1e3.
        shots : int
            The number of shots per point, default is 100.
        sweep : float
            The full sweep (span) width in GHz, default is 1 GHz. 
            Note that the sweep is centered around the freq frequency, so the sweep goes from freq-sweep/2 to freq+sweep/2.
        step: float
            The frequency step in GHz, default 0.05
        dim: int
            The number of points in the sweep, default 25.

        

        """

        name = 'LO Sweep Tune Sequence'
        B = kwargs.pop('B', None)
        reptime = kwargs.pop('reptime', 1e3) # 0.1ms
        shots=kwargs.pop('shots', 100)
        averages=kwargs.pop('averages', 1)

        if all(key in kwargs for key in ('sweep', 'step', 'dim')):
            raise ValueError("Cannot specify sweep, step and dim together. Specify only two of them.")
        elif 'sweep' not in kwargs:
            step = kwargs.pop('step', 0.05)
            dim = int(kwargs.pop('dim', 25))
            sweep = step*(dim-1)
        elif 'step' not in kwargs:
            sweep = kwargs.pop('sweep', 1)
            dim = int(kwargs.pop('dim', 25))
            step = sweep/(dim-1)
        elif 'dim' not in kwargs:
            sweep = kwargs.pop('sweep', 1)
            step = kwargs.pop('step', 0.05)
            dim = int(sweep/step)+1
        

        super().__init__(name=name, B=B, freq=freq, reptime=reptime,shots=shots,averages=averages, **kwargs)

        self.freq = Parameter(name='freq', value=freq-sweep/2, step=step, dim=dim, unit='GHz', virtual=False)
        self._build()
        self.evolution([self.freq])

    def _build(self):
        self.addPulse(RectPulse(t=350,tp=300, freq=0, scale=1, flipangle='hard'))
        self.addPulse(Detection(t=512,tp=1024))
