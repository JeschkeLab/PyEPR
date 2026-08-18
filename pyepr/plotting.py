import matplotlib.pyplot as plt
import numpy as np
import pyepr as epr

def plot_sequence_freq(sequence, resonator=None, spectrum=None, title=None, fig=None, axs=None, **kwargs):
    """
    Plots the relative pulses in the frequency domain. If a resonator or spectrum is provided, these are plotted as well.  

    Parameters
    ----------
    sequence : epr.Sequence
        The pulse sequence to be plotted.
    resonator : epr.ResonatorProfile, optional
        The resonator profile to be plotted. If provided, the resonator profile will be plotted on a secondary y-axis.
    spectrum : epr.FieldSweepAnalysis, optional
        The field sweep spectrum to be plotted. If provided, the spectrum will be plotted on the same axis as the pulse sequence.
    title : str, optional
        The title of the plot. If not provided, no title will be set.
    fig : matplotlib.figure.Figure, optional
        The figure to plot on. If not provided, a new figure will be created.
    axs : matplotlib.axes.Axes, optional
        The axes to plot on. If not provided, new axes will be created.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot.
    """

    if fig is None or axs is None:
        figsize = kwargs.get('figsize', (6,4))
        fig, axs = plt.subplots(figsize=figsize, layout='constrained')
    elif fig is None and isinstance(axs,plt.Axes):
        fig = axs.figure

    fc = sequence.freq.value
    B = sequence.B.value
    if resonator is not None:
        axs_res = axs.twinx()
        axs_res.plot(resonator.freqs, resonator.profile*1e3,'.-', label='Resonator Profile', color=epr.primary_colors[0])
        axs_res.set_ylabel(r'Resonator Profile $\nu_1$ / MHz', color=epr.primary_colors[0])
        BW = resonator.fc/resonator.q
        freqs = np.linspace(-BW*3, BW*3, 100)
    else:
        axs_res = None
        freqs = np.linspace(-0.5, 0.5, 100)

    if spectrum is not None and hasattr(spectrum, 'gyro'):
        gyro = spectrum.gyro
        freq_peak = B*gyro
        axs.plot(spectrum.fs_x+freq_peak,spectrum.data, color=epr.primary_colors[1], label='Field Sweep Spectrum')

    for i,pulse in enumerate(sequence.pulses):
        if isinstance(pulse,epr.Detection):
            continue
        pulse:epr.Pulse

        profile = pulse.exciteprofile(freqs=freqs)[:,2]
        rescaled_profile = (profile *-0.5 + 0.5) / (pulse.flipangle.value/np.pi)
        label = f"Pulse {i}"

        axs.fill_between(freqs+fc, rescaled_profile, label=label,alpha=0.3,color=epr.secondary_colors[i%len(epr.secondary_colors)])

    # Cretate parallel x axis offset by seq.freq.value
    secax = axs.secondary_xaxis('top', functions=(lambda x: x - fc, lambda x: x + fc))
    secax.set_xlabel('Pulse Frequency / GHz')  # adjust label/units as needed
    axs.set_xlabel('Lab Frequency / GHz')
    axs.set_ylabel('Normalized Pulse Efficiency ')
    axs.set_xlim(freqs[0]+fc, freqs[-1]+fc)

    lines1, labels1 = axs.get_legend_handles_labels()
    if axs_res is not None:
        lines2, labels2 = axs_res.get_legend_handles_labels()
        lines = lines2+ lines1 
        labels = labels2+ labels1
    else:
        lines = lines1
        labels = labels1
    fig.legend(lines, labels, ncol=3, loc='outside lower left', fontsize=8)

    if title is not None:
        fig.suptitle(title)

    return fig

def plot_sequence_time(sequence, title=None, fig=None, axs=None):
    """
    Creates a graphical summary of a pulse sequence over time. 
    Showing the pulses at their respective times and durations, as well as the sequence of events in the experiment.
    
    """

    pass