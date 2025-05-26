from matplotlib.figure import Figure
import matplotlib.cm as cm
import numpy as np
from deerlab import noiselevel
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pyepr.sequences import Sequence
import deerlab as dl
from scipy.linalg import svd
from pyepr.colors import primary_colors
# ===========================================================================


class CarrPurcellAnalysis():

    def __init__(self, dataset, sequence: Sequence = None) -> None:
        """Analysis and calculation of Carr Purcell decay. 

        Parameters
        ----------
        dataset : 
            _description_

        Attributes
        ----------
        axis : xr.DataArray
            The time axis representing the interpulse delay.
        """
        # self.axis = dataset.axes[0]
        # self.data = dataset.data
        if 'tau1' in dataset.coords:
            self.axis = dataset['tau1']
        elif 'tau' in dataset.coords:
            self.axis = dataset['tau']
        elif 't' in dataset.coords:
            self.axis = dataset['t']
        elif 'step' in dataset.coords:
            self.axis = dataset['step'] 
        else:
            self.axis = dataset['X']
        
        dataset = dataset.epr.correctphasefull
        self.data = dataset.data.real
        self.dataset = dataset

        data = self.data / np.max(self.data)
        self.noise = noiselevel(data)
        # if sequence is None and hasattr(dataset,'sequence'):
        #     self.seq = dataset.sequence
        # else:
        #     self.seq = sequence
        pass
    
    def fit(self, type: str = "mono",**kwargs):
        """Fit the experimental CP decay

        Parameters
        ----------
        type : str, optional
            Either a mono or double exponential decay model, by default "mono"

        """

        data = self.data
        data /= np.max(data)

        # if type == "mono":
        #     self.func = lambda x, a, b, e: a*np.exp(-b*x**e)
        #     p0 = [1, 1, 2]
        #     bounds = ([0, 0, 0],[2, 1000, 10])
        # elif type == "double":
        #     self.func = lambda x, a, b, e, c, d, f: a*np.exp(-b*x**e) + c*np.exp(-d*x**f)
        #     p0 = [1, 1, 2, 1, 1, 2]
        #     bounds = ([0, 0, 0, 0, 1, 0],[2, 1000, 10, 2, 1000, 10])
        # else:
        #     raise ValueError("Type must be one of: mono")s
        
        # self.fit_type = type
        # self.fit_result = curve_fit(self.func, self.axis, data, p0=p0, bounds=bounds)
        monoModel = dl.bg_strexp
        monoModel.name = 'Stretched exponential'
        doubleModel = dl.bg_sumstrexp
        doubleModel.weight1.ub = 500
        doubleModel.name = "Sum of two stretched exponentials"

        testModels = []
        if type == "mono":
            testModels.append(monoModel)
        elif type == "double":
            testModels.append(doubleModel)

        else: # type == "auto"
            testModels = [monoModel, doubleModel]

        results = []
        for model in testModels:
            tmp_result = dl.fit(model,data,self.axis,reg=False,**kwargs)
            results.append(tmp_result)
        
        if len(results) == 1:
            self.fit_result = results[0]
            self.fit_model = testModels[0]
        else:
            # Select based of R2 and weight of first component
            R2 = [result.stats['R2'] for result in results]
            self.fit_result = results[np.argmax(R2)]
            self.fit_model = testModels[np.argmax(R2)]
            print(f"Selected model: {self.fit_model.description}")
        
        return self.fit_result

    def plot(self, norm: bool = True, ci=50, axs=None, fig=None) -> Figure:
        """Plot the carr purcell decay with fit, if avaliable.

        Parameters
        ----------
        norm : bool, optional
            Normalise the fit to a maximum of 1, by default True
        ci : int, optional
            The percentage confidence interval to plot, by default 50
        

        Returns
        -------
        Figure
            The figure.
        """

        if norm is True:
            data = self.data
            data /= np.max(data)

        if axs is None and fig is None:
            fig, axs = plt.subplots()

        if hasattr(self, "fit_result"):
            x = self.axis
            V = self.fit_result.evaluate(self.fit_model, x)*self.fit_result.scale
            fitUncert = self.fit_result.propagate(self.fit_model, x)
            VCi = fitUncert.ci(ci)*self.fit_result.scale
            ub = VCi[:,1]
            lb = VCi[:,0]
            # ub = self.fit_model(x,*self.fit_result.paramUncert.ci(ci)[:-1,1])*self.fit_result.paramUncert.ci(ci)[-1,1]
            # lb = self.fit_model(x,*self.fit_result.paramUncert.ci(ci)[:-1,0])*self.fit_result.paramUncert.ci(ci)[-1,0]
            axs.plot(self.axis, data, '.', label='data', color='0.6', ms=6)
            axs.plot(x, V, label='fit', color=primary_colors[0], lw=2)
            if ci is not None:
                axs.fill_between(x, lb, ub, color=primary_colors[0], alpha=0.3, label=f"{ci}% CI")

            axs.legend()
        else:
            axs.plot(self.axis, data, label='data')

        axs.set_xlabel('Time / us')
        axs.set_ylabel('Normalised Amplitude')
        return fig
    
    def check_decay(self,level=0.1):
        """
        Checks that the data has decayed by over 90% in the first half, and less than 90% in the first quarter.

        Parameters
        ----------
        level : float, optional
            The level to check the decay, by default 0.05

        Returns
        -------
        int
            0 if both conditions are met, 1 if a longer decay is needed, and -1 if the decay is too long.
        
        """
        n_points = len(self.axis)
        if hasattr(self,"fit_result"):
            # decay = self.func(self.axis, *self.fit_result[0]).data
            x = self.axis
            decay = self.fit_result.evaluate(self.fit_model, x)*self.fit_result.scale
            if (decay[:int(n_points*0.50)].min() < level) and (decay[:int(n_points*0.25)].min() > level):
                return 0
            elif decay[:int(n_points*0.25)].min() < level:
                return 1
            elif decay[:int(n_points*0.50)].min() > level:
                return -1
        else:
            raise ValueError("No fit result found")

    def find_optimal(
            self, SNR_target, target_time: float, target_step, averages=None) -> float:
        """Calculate the optimal inter pulse delay for a given total measurment
        time. 

        Parameters
        ----------
        SNR_target: float,
            The Signal to Noise ratio target.
        target_time : float
            The target time in hours
        target_shrt : float
            The shot repettition time of target in seconds
        target_step: float
            The target step size in ns.
        averages : int, optional
            The total number of shots taken, by default None. If None, the
            number of shots will be calculated from the dataset.
        

        Returns
        -------
        float
            The calculated optimal time in us
        """
        # time_per_point = shrt * averages
        dataset = self.dataset
        if averages is None:
            averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
        target_shrt = dataset.reptime * 1e-6

        data = np.abs(self.data)
        data /= np.max(data)

        if hasattr(self,"fit_result"):
            calc_data = self.func(self.axis.data,*self.fit_result[0])
        else:
            calc_data = data

        # averages = self.seq.shots.value * self.seq.averages.value
        self.noise = noiselevel(data)
        data_snr = calc_data / self.noise
        data_snr_avgs = data_snr / np.sqrt(averages)

        # Target time
        target_time = target_time * 3600
        target_step_us = target_step * 1e-3
        g = (target_time * target_step / target_shrt) * 1/(self.axis.data)
        f = (SNR_target/data_snr_avgs)**2

        self.optimal = self.axis.data[np.argmin(np.abs(g-f))]
        return self.optimal
    
    def __call__(self, x, norm=True, SNR=False, source=None):
        """
        Evaluate the fit or data at a given x value.

        Parameters
        ----------
        x : float
            The x value to evaluate the data at.
        norm : bool, optional
            Normalise the data to the maximum, by default True
        SNR : bool, optional
            Return the SNR_per_sqrt(shot) for this data point, by default False
            If True, the data is normalised to the maximum of the data.
        source : str, optional
            The source of the data, either 'fit' or 'data', by default None
            If None, the source is determined by the presence of a fit result.
        
        """
        
        if source is 'fit' or (source is None and hasattr(self,'fit_result')):
            V = self.fit_result.evaluate(self.fit_model, x)*self.fit_result.scale
            if not norm and SNR is not True: # Fit is normalised to 1 by default
                V *= np.max(self.data)
        elif source is 'data' or (source is None and not hasattr(self,'fit_result')):
            x_idx = np.abs(self.axis - x).argmin()
            V = self.data[x_idx]

            if norm is True or SNR is True:
                V /= np.max(self.data)

        if SNR is True:
            V /= self.noise
            V /= np.sqrt(self.dataset.nAvgs * self.dataset.shots * self.dataset.nPcyc)

        # return single value if x is a single value
        if np.isscalar(x):
            return V[0]
        else:
            return V

class HahnEchoRelaxationAnalysis():

    def __init__(self, dataset) -> None:
        """Analysis, fitting and plotting for the HahnEchoRelaxation Sequence. 

        Parameters
        ----------
        dataset : xarray.DataArray
            The dataset to be analysed, with the time axis contained.

        Attributes
        ----------
        axis : xr.DataArray
            The time axis representing the interpulse delay.
        """
        # self.axis = dataset.axes[0]
        # self.data = dataset.data
        if 'tau1' in dataset.coords:
            self.axis = dataset['tau1']
        elif 'tau' in dataset.coords:
            self.axis = dataset['tau']
        elif 't' in dataset.coords:
            self.axis = dataset['t']
        elif 'step' in dataset.coords:
            self.axis = dataset['step'] 
        else:
            self.axis = dataset['X']
        
        dataset = dataset.epr.correctphasefull
        self.data = dataset.data.real
        self.dataset = dataset

        data = self.data / np.max(self.data)
        self.noise = noiselevel(data)

        
        pass
    
    def fit(self, type: str = "mono",**kwargs):
        """Fit the experimental CP decay

        Parameters
        ----------
        type : str, optional
            Either a mono or double exponential decay model, by default "mono"

        """

        data = self.data
        data /= np.max(data)
        monoModel = dl.bg_strexp
        monoModel.name = 'Stretched exponential'
        doubleModel = dl.bg_sumstrexp
        doubleModel.weight1.ub = 200
        doubleModel.name = "Sum of two stretched exponentials"

        testModels = []
        if type == "mono":
            testModels.append(monoModel)
        elif type == "double":
            testModels.append(doubleModel)

        else: # type == "auto"
            testModels = [monoModel, doubleModel]

        results = []
        for model in testModels:
            results.append(dl.fit(model,data,self.axis,reg=False,**kwargs))
        
        if len(results) == 1:
            self.fit_result = results[0]
            self.fit_model = testModels[0]
        else:
            # Select based of R2
            R2 = [result.stats['R2'] for result in results]
            self.fit_result = results[np.argmax(R2)]
            self.fit_model = testModels[np.argmax(R2)]
            print(f"Selected model: {self.fit_model.description}")
        
        return self.fit_result

    def plot(self, norm: bool = True, ci=50, axs=None, fig=None) -> Figure:
        """Plot the carr purcell decay with fit, if avaliable.

        Parameters
        ----------
        norm : bool, optional
            Normalise the fit to a maximum of 1, by default True
        ci : int, optional
            The percentage confidence interval to plot, by default 50
        

        Returns
        -------
        Figure
            The figure.
        """

        if norm is True:
            data = self.data
            data /= np.max(data)

        if axs is None and fig is None:
            fig, axs = plt.subplots()

        if hasattr(self, "fit_result"):
            x = self.axis
            V = self.fit_result.evaluate(self.fit_model, x)*self.fit_result.scale
            fitUncert = self.fit_result.propagate(self.fit_model, x)
            VCi = fitUncert.ci(ci)*self.fit_result.scale
            ub = VCi[:,1]
            lb = VCi[:,0]
            # ub = self.fit_model(x,*self.fit_result.paramUncert.ci(ci)[:-1,1])*self.fit_result.paramUncert.ci(ci)[-1,1]
            # lb = self.fit_model(x,*self.fit_result.paramUncert.ci(ci)[:-1,0])*self.fit_result.paramUncert.ci(ci)[-1,0]
            axs.plot(self.axis, data, '.', label='data', color='0.6', ms=6)
            axs.plot(x, V, label='fit', color=primary_colors[0], lw=2)
            if ci is not None:
                axs.fill_between(x, lb, ub, color=primary_colors[0], alpha=0.3, label=f"{ci}% CI")

            axs.legend()
        else:
            axs.plot(self.axis, data, label='data')

        axs.set_xlabel('Time / us')
        axs.set_ylabel('Normalised Amplitude')
        return fig
    
    def check_decay(self,level=0.1):
        """
        Checks that the data has decayed by over 90% in the first half, and less than 90% in the first quarter.

        Parameters
        ----------
        level : float, optional
            The level to check the decay, by default 0.1

        Returns
        -------
        int
            0 if both conditions are met, 1 if a longer decay is needed, and -1 if the decay is too long.
        
        """
        n_points = len(self.axis)
        if hasattr(self,"fit_result"):
            # decay = self.func(self.axis, *self.fit_result[0]).data
            x = self.axis
            decay = self.fit_result.evaluate(self.fit_model, x)*self.fit_result.scale
            if (decay[:int(n_points*0.50)].min() < level) and (decay[:int(n_points*0.25)].min() > level):
                return 0
            elif decay[:int(n_points*0.25)].min() < level:
                return 1
            elif decay[:int(n_points*0.50)].min() > level:
                return -1
        else:
            raise ValueError("No fit result found")

    def __call__(self, x, norm=True, SNR=False, source=None):
        """
        Evaluate the fit or data at a given x value.

        Parameters
        ----------
        x : float
            The x value to evaluate the data at.
        norm : bool, optional
            Normalise the data to the maximum, by default True
        SNR : bool, optional
            Return the SNR_per_sqrt(shot) for this data point, by default False
        source : str, optional
            The source of the data, either 'fit' or 'data', by default None
            If None, the source is determined by the presence of a fit result.
        
        """
        
        if source is 'fit' or (source is None and hasattr(self,'fit_result')):
            V = self.fit_result.evaluate(self.fit_model, x)*self.fit_result.scale
            if not norm: # Fit is normalised to 1 by default
                V *= np.max(self.data)
        elif source is 'data' or (source is None and not hasattr(self,'fit_result')):
            x_idx = np.abs(self.axis - x).argmin()
            V = self.data[x_idx]

            if norm is True:
                V /= np.max(self.data)

        if SNR is True:
            V /= self.noise
            V /= np.sqrt(self.dataset.nAvgs * self.dataset.shots * self.dataset.nPcyc)

        # return single value if x is a single value
        if np.isscalar(x):
            return V[0]
        else:
            return V
        
class ReptimeAnalysis():

    def __init__(self, dataset, sequence: Sequence = None) -> None:
        """Analysis and calculation of Reptime based saturation recovery. 

        Parameters
        ----------
        dataset :
            The dataset to be analyzed.
        sequence : Sequence, optional
            The sequence object describing the experiment. (not currently used)
        """
        # self.axis = dataset.axes[0]
        self.axis = dataset['reptime']
        # if self.axis.max() > 1e4:
        #     self.axis /= 1e3 # ns -> us
        # self.data = dataset.data/np.max(dataset.data)
        
        if np.iscomplexobj(dataset.data):
            self.data = dataset.epr.correctphase
        else:
            self.data = dataset

        self.data.data /= np.max(self.data.data)
        self.seq = sequence
        pass

    def fit(self,type='SE', **kwargs):

        if type == 'SE': # stetch exponential recovery
            def func(t,A,T1,xi):
                return A*(1-np.exp(-(t/T1)**xi))
            p0 = [1,1.8e3,1]
        elif type.lower() == 'exp': # exponential recovery
            def func(t,A,T1):
                return A*(1-np.exp(-t/T1))
            p0 = [1,1.8e3]
        self.func = func

        if 'p0' in kwargs:
            p0 = kwargs.pop('p0')
        # mymodel = dl.Model(func,constants='t')
        # mymodel.T1.set(lb=0,ub=np.inf,par0=1.8e3)
        # mymodel.T1.unit = 'us'
        # mymodel.T1.description = 'T1 time'

        # results = dl.fit(mymodel,self.data.real,self.axis,reg=False,**kwargs)
        # self.fit_result = results

        self.fit_result = curve_fit(func, self.axis, self.data, p0=p0,**kwargs)

        return self.fit_result

    def plot(self, axs=None, fig=None):

        if axs is None and fig is None:
            fig, axs = plt.subplots()

        # if hasattr(self,'fit_result'):
        #     # renormalise data to fit amplitude
        #     data = self.data/self.fit_result[0][0]
        # else:
        data = self.data

        axs.plot(self.axis/1e3, data, '.', label='data', color='0.6', ms=6)
        
        if hasattr(self,'fit_result'):
            axs.plot(self.axis/1e3, self.func(self.axis,*self.fit_result[0]), label='Fit', color=primary_colors[0], lw=2)
            axs.set_xlim(*axs.get_xlim())
            axs.set_ylim(*axs.get_ylim())
            ylim = axs.get_ylim()
            axs.vlines(self.fit_result[0][1]/1e3,*ylim,linestyles='dashed',label='T1 = {:.3g} ms'.format(self.fit_result[0][1]/1e3),colors=primary_colors[1])

            if hasattr(self,'optimal'):
                axs.vlines(self.optimal/1e3,*ylim,linestyles='dashed',label='Optimal = {:.3g} ms'.format(self.optimal/1e3),colors=primary_colors[2])

        axs.set_xlabel('Reptime / ms')
        axs.set_ylabel('Normalised signal')
        axs.legend()
        return fig

    def calc_optimal_reptime(self, recovery=0.9):
        # Calculates the x% recovery time
        if recovery is not None:
            T1 = self.fit_result[0][1]
            if self.fit_result[0].shape[0] == 3:
                xi = self.fit_result[0][2]
            else:
                xi = 1
            self.optimal = T1 * np.log(1/(1-recovery))**(1/xi)            
        else:
            t = self.axis
            optimal_vals = self.func(t,*self.fit_result[0])* 1/np.sqrt(t)
            self.optimal = t[np.nanargmax(optimal_vals)]
        return self.optimal

def detect_ESEEM(dataset,type='deuteron', threshold=1.5):
    """Detect if the dataset is an ESEEM experiment.

    Parameters
    ----------
    dataset : xr.DataArray
        The dataset to be analyzed.
    
    type : str, optional
        The type of ESEEM experiment, either deuteron or proton, by default 'deuteron'
    
    threshold : float, optional
        The SNR threshold for detection, by default 1.5

    Returns
    -------
    bool
        True if ESEEM is detected, False if not.
    """
    

    D_freq = 4.10663 * dataset.B *1e-4 *np.pi /2
    P_freq = 26.75221 * dataset.B *1e-4 *np.pi /2

    def find_pnl(freq):
        fft_data = np.abs(dataset.epr.fft)
        index = np.abs(fft_data.X - freq).argmin().data

        peak = 2 /fft_data.size * fft_data[index]

        noiselevel = 2/fft_data.size * fft_data[index-8:index+8].mean()

        return peak/noiselevel
        
    if type == 'deuteron':
        peak = find_pnl(D_freq)
    elif type == 'proton':
        peak = find_pnl(P_freq)
    else:
        raise ValueError('type must be deuteron or proton')

    if peak > threshold:
        return True
    else:
        return False

cmap = ['#D95B6F','#42A399']

def plot_1Drelax(*args,fig=None, axs=None,cmap=cmap):
    """
    Create a superimposed plot of relaxation data and fits.

    Parameters
    ----------
    args : ad.Analysis
        The 1D relaxation data to be plotted.

    fig : Figure, optional
        The figure to plot to, by default None
    axs : Axes, optional
        The axes to plot to, by default None
    cmap : list, optional
        The color map to use, by default ad.cmap
    
    """

    if fig is None and axs is None:
        fig, axs = plt.subplots(1,1, figsize=(5,5))
    elif axs is None:
        axs = fig.subplots(1,1)

    for i,arg in enumerate(args): 
        if arg.dataset.seq_name == 'T2RelaxationSequence':
            xscale = 2
            label='Hahn Echo'
        elif arg.dataset.seq_name == 'RefocusedEcho1DSequence':
            xscale = 2
            label='1D Refocused Echo'
        elif arg.dataset.seq_name == 'CarrPurcellSequence':
            xscale = 4
            label='CP-2'
        elif (arg.dataset.seq_name == 'DEERSequence') or (arg.dataset.seq_name == '5pDEER'):
            xscale = 4
            label='CP-2'

        else:
            xscale = 4
            label='CP-2'
        
        axs.plot(arg.axis*xscale, arg.data/arg.data.max(), '.', label=label,alpha=0.5,color=cmap[i],mec='none')
        if hasattr(arg, 'func'):
            print('The scipy fitting elements are being deprecated, please use DeerLab fitting')
            V = arg.func(arg.axis,*arg.fit_result[0])
        else:
            V = arg.fit_model(arg.axis,*arg.fit_result.param[:-1])*arg.fit_result.scale
        axs.plot(arg.axis*xscale, V, '-',alpha=1,color=cmap[i], lw=2)

    axs.legend()
    axs.set_xlabel('Total Sequence Length / $\mu s$')
    axs.set_ylabel('Signal / $ A.U. $')

    return fig


            

        