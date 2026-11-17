import matplotlib.pyplot as plt
import numpy as np
import pyepr as epr
import numbers

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

def draw_sequence(sequence, ax=None,**kwargs):
        """
        Draws the sequence as an anotated matplotlib figure, showing 
        pictorial representations of the pulses and delays in the sequence. It
        does not show the phase cycle and nothing is to scale.

        Pulses are drawn with their amplitude envelope and labelled with their
        flip angle. Frequency swept pulses are additionally shown with a chirp
        of increasing (or decreasing) frequency. Pulses at a frequency other
        than the detection frequency are drawn on a separate channel below.
        The delay between each consecutive event is annotated beneath the
        sequence.

        Evolution axes are colour coded and listed beneath the sequence.
        Delays and pulse parameters that change along an evolution axis are
        drawn in that axis colour, and labelled with the range they cover.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axes to draw onto. If None, a new figure is created.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure containing the drawing.
        """
        from pyepr.colors import primary_colors, secondary_colors

        time_scales = {"ns": 1, "us": 1e3, "ms": 1e6, "s": 1e9}

        def to_ns(param):
            return param.value * time_scales.get(param.unit, 1)

        evo_params = getattr(sequence, "evo_params", None)
        if evo_params is None:
            raise ValueError("The sequence has no evolution parameters to draw. Please call Sequence.evolution first.")
        evo_uuids = [param.uuid for param in evo_params]
        evo_colors = {uid: secondary_colors[i % len(secondary_colors)]
                      for i, uid in enumerate(evo_uuids)}

        def evo_axes(param):
            # Offsets from the initial value, for each evolution axis the
            # parameter follows
            axes = {}
            for axis in getattr(param, "axis", []):
                if axis["uuid"] in evo_uuids:
                    axes[axis["uuid"]] = np.asarray(axis["axis"], dtype=float)
            return axes

        def fmt_value(value, unit):
            if unit in time_scales:
                return fmt_delay(value * time_scales[unit])
            if unit in (None, "", "None"):
                return f"{value:.4g}"
            return f"{value:.5g} {unit}"

        def fmt_range(start, end, unit):
            return f"{fmt_value(start, unit)} → {fmt_value(end, unit)}"

        def pulse_evolution(pulse):
            # Evolving pulse parameters, other than the start time
            lines = []
            for name, param in vars(pulse).items():
                if name == "t" or not isinstance(param, epr.Parameter) or param.value is None:
                    continue
                axes = evo_axes(param)
                if len(axes) == 0:
                    continue
                start = param.value + sum(axis[0] for axis in axes.values())
                end = param.value + sum(axis[-1] for axis in axes.values())
                lines.append((f"{name}: {fmt_range(start, end, param.unit)}",
                              evo_colors[next(iter(axes))]))
            return lines

        def centre_freq(pulse):
            if hasattr(pulse, "freq"):
                return pulse.freq.value
            elif hasattr(pulse, "init_freq") and hasattr(pulse, "final_freq"):
                return (pulse.init_freq.value + pulse.final_freq.value) / 2
            return 0

        def is_swept(pulse):
            return hasattr(pulse, "init_freq") and hasattr(pulse, "final_freq")

        def flip_label(pulse):
            flipangle = getattr(pulse, "flipangle", None)
            if flipangle is None or flipangle.value is None:
                return ""
            if isinstance(flipangle.value, str):
                return flipangle.value
            frac = flipangle.value / np.pi
            for den in range(1, 9):
                num = int(np.round(frac * den))
                if num != 0 and np.isclose(num / den, frac):
                    num_str = {1: "", -1: "-"}.get(num, str(num))
                    den_str = "" if den == 1 else f"/{den}"
                    return rf"${num_str}\pi{den_str}$"
            return f"{np.degrees(flipangle.value):.0f}°"

        def pulse_width(pulse):
            # Pictorial widths, a pi pulse is twice as wide as a pi/2 pulse
            if isinstance(pulse, epr.Detection) or is_swept(pulse):
                return 1.2
            flipangle = getattr(pulse, "flipangle", None)
            if flipangle is None or not isinstance(flipangle.value, numbers.Number):
                width = 0.5
            else:
                width = float(np.clip(0.7 * abs(flipangle.value) / np.pi, 0.25, 1.0))
            if isinstance(pulse, (epr.GaussianPulse, epr.SincPulse)):
                width = max(width, 0.8)  # Shaped pulses need room to show their envelope
            return width

        def envelope(pulse):
            # Idealised amplitude envelope over x in [-1, 1]
            x = np.linspace(-1, 1, 400)
            if isinstance(pulse, epr.GaussianPulse):
                AM = np.exp(-(x / 0.45)**2)
            elif isinstance(pulse, epr.SincPulse):
                AM = np.sinc(3 * x)
            elif isinstance(pulse, epr.HSPulse):
                AM = 1 / np.cosh(4 * np.abs(x)**pulse.order1.value)
            elif isinstance(pulse, epr.WURSTPulse):
                AM = 1 - np.abs(np.sin(np.pi * x / 2))**pulse.order.value
            elif isinstance(pulse, epr.ChirpPulse):
                AM = np.clip((1 - np.abs(x)) / 0.15, 0, 1)
                AM = np.sin(np.pi / 2 * AM)**2
            else:
                return np.array([-1, -1, 1, 1]), np.array([0, 1, 1, 0])
            AM[[0, -1]] = 0
            return x, AM

        def fmt_delay(dt):
            if abs(dt) >= 1e3:
                return f"{dt/1e3:.4g} µs"
            return f"{dt:.4g} ns"

        # Collect the events in time order
        events = [pulse for pulse in sequence.pulses
                  if pulse.t is not None and pulse.t.value is not None]
        if len(events) == 0:
            raise ValueError("The sequence contains no pulses to draw.")
        events.sort(key=lambda pulse: to_ns(pulse.t))

        # Assign each event to a frequency channel, the observer channel is
        # the detection frequency
        dets = [pulse for pulse in events if isinstance(pulse, epr.Detection)]
        obs_freq = centre_freq(dets[0] if len(dets) > 0 else events[0])
        channel_freqs = [obs_freq]
        channels = []
        for pulse in events:
            freq = centre_freq(pulse)
            match = [i for i, f in enumerate(channel_freqs) if np.isclose(f, freq, atol=1e-4)]
            if len(match) == 0:
                channel_freqs.append(freq)
                match = [len(channel_freqs) - 1]
            channels.append(match[0])

        # Group simultaneous events into slots, evenly spaced along the x-axis
        slot_times = []
        slot_axes = []
        slot_of_event = []
        for pulse in events:
            t = to_ns(pulse.t)
            if len(slot_times) == 0 or not np.isclose(t, slot_times[-1], atol=1e-3):
                slot_times.append(t)
                scale = time_scales.get(pulse.t.unit, 1)
                slot_axes.append({uid: axis * scale for uid, axis in evo_axes(pulse.t).items()})
            slot_of_event.append(len(slot_times) - 1)

        def delay_evolution(i):
            # Change of the delay between slot i and i+1 along each evolution axis
            diffs = {}
            for uid in evo_uuids:
                diff = slot_axes[i+1].get(uid, 0) - slot_axes[i].get(uid, 0)
                if np.ndim(diff) > 0 and not np.allclose(diff, 0):
                    diffs[uid] = diff
            return diffs

        slot_widths = np.zeros(len(slot_times))
        for pulse, slot in zip(events, slot_of_event):
            slot_widths[slot] = max(slot_widths[slot], pulse_width(pulse))

        slot_x = np.zeros(len(slot_times))
        for i in range(1, len(slot_times)):
            slot_x[i] = slot_x[i-1] + slot_widths[i-1]/2 + slot_widths[i]/2 + 1.8

        # Layout
        height = 1.0
        channel_spacing = 2.2
        baselines = [-i * channel_spacing for i in range(len(channel_freqs))]
        # Leave room under the bottom channel for evolving pulse parameters
        n_bottom_lines = max([len(pulse_evolution(pulse)) for pulse, channel
                              in zip(events, channels) if channel == len(channel_freqs) - 1])
        arrow_y = baselines[-1] - 0.6 - 0.3 * n_bottom_lines
        legend_y = arrow_y - 0.75
        y_min = legend_y - 0.35 * len(evo_params) + (0.2 if len(evo_params) else 0.15)
        x_min = slot_x[0] - slot_widths[0]/2 - 1.0
        x_max = slot_x[-1] + slot_widths[-1]/2 + 1.0

        if ax is None:
            fig_width = float(np.clip(0.75 * (x_max - x_min) + 1.5, 6, 18))
            fig_height = 0.7 * (baselines[0] + height + 0.6 - y_min) + 0.6
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        else:
            fig = ax.figure

        colors = [primary_colors[i % len(primary_colors)] for i in range(len(channel_freqs))]

        # Channel baselines and labels
        for i, (y0, freq) in enumerate(zip(baselines, channel_freqs)):
            ax.plot([x_min, x_max], [y0, y0], color='k', lw=1.2, zorder=1)
            if len(channel_freqs) == 1:
                label = "mw"
            elif i == 0:
                label = r"$\nu_{obs}$"
            elif len(channel_freqs) == 2:
                label = r"$\nu_{pump}$"
            else:
                label = rf"$\nu_{{{i+1}}}$"
            ax.text(x_min - 0.2, y0 + 0.3*height, label, ha='right', va='center', fontsize=12)
            if i > 0:
                offset = (freq - obs_freq) * 1e3
                ax.text(x_min - 0.2, y0 - 0.15*height, rf"$\Delta\nu$ = {offset:+.0f} MHz",
                        ha='right', va='center', fontsize=8, color='0.35')

        # Vertical guides from each event down to the delay arrows
        for x in slot_x:
            ax.plot([x, x], [arrow_y, baselines[0]], color='0.75', lw=0.8, ls=':', zorder=0)

        # Pulses and detection events
        for pulse, slot, channel in zip(events, slot_of_event, channels):
            x0 = slot_x[slot]
            y0 = baselines[channel]
            w = pulse_width(pulse)
            color = colors[channel]

            for j, (text, evo_color) in enumerate(pulse_evolution(pulse)):
                ax.text(x0, y0 - 0.12 - 0.3*j, text, ha='center', va='top',
                        fontsize=8, color=evo_color)

            if isinstance(pulse, epr.Detection):
                ax.add_patch(plt.Rectangle(
                    (x0 - w/2, y0), w, height, facecolor=color, alpha=0.1,
                    edgecolor=color, ls='--', lw=0.8, zorder=1))
                xs = np.linspace(-1, 1, 400)
                env = np.exp(-(xs / 0.35)**2)
                ax.fill_between(x0 + xs*w/2, y0, y0 + 0.8*height*env,
                                color=color, alpha=0.25, lw=0, zorder=2)
                ax.plot(x0 + xs*w/2, y0 + 0.8*height*env*np.cos(2*np.pi*2.5*xs),
                        color=color, lw=1.2, zorder=3)
                label = "echo"
                if hasattr(pulse, "n_echoes"):
                    label += rf" $\times${pulse.n_echoes.value}"
                ax.text(x0, y0 + height + 0.1, label, ha='center', va='bottom', fontsize=10)
                continue

            xs, AM = envelope(pulse)
            px = x0 + xs*w/2
            ax.fill_between(px, y0, y0 + height*AM, color=color, alpha=0.85 if not is_swept(pulse) else 0.35,
                            lw=0, zorder=2)
            ax.plot(px, y0 + height*AM, color=color, lw=1.0, zorder=3)

            if is_swept(pulse):
                # Draw a chirp inside the envelope to show the frequency sweep
                xs_fine = np.linspace(-1, 1, 600)
                AM_fine = np.interp(xs_fine, xs, AM)
                f0, f1 = (2, 9) if pulse.final_freq.value >= pulse.init_freq.value else (9, 2)
                u = (xs_fine + 1) / 2
                phase = 2*np.pi*(f0*u + (f1 - f0)*u**2/2)
                ax.plot(x0 + xs_fine*w/2,
                        y0 + 0.5*height*AM_fine*(1 + 0.85*np.sin(phase)),
                        color=color, lw=0.8, zorder=3)

            ax.text(x0, y0 + height*max(AM.max(), 0) + 0.1, flip_label(pulse),
                    ha='center', va='bottom', fontsize=12)

        # Delays between consecutive events, evolving delays are coloured by
        # their evolution axis and given as a range
        for i in range(len(slot_times) - 1):
            dt = slot_times[i+1] - slot_times[i]
            diffs = delay_evolution(i)
            if len(diffs) > 0:
                color = evo_colors[next(iter(diffs))]
                label = fmt_range(dt + sum(diff[0] for diff in diffs.values()),
                                  dt + sum(diff[-1] for diff in diffs.values()), "ns")
                lw = 1.6
            else:
                color, label, lw = '0.3', fmt_delay(dt), 0.9
            ax.annotate("", xy=(slot_x[i+1], arrow_y), xytext=(slot_x[i], arrow_y),
                        arrowprops=dict(arrowstyle='<->', color=color, lw=lw,
                                        shrinkA=0, shrinkB=0))
            ax.text((slot_x[i] + slot_x[i+1]) / 2, arrow_y - 0.12, label,
                    ha='center', va='top', fontsize=9,
                    color='k' if len(diffs) == 0 else color)

        # Legend of the evolution axes
        for i, param in enumerate(evo_params):
            axis = evo_axes(param).get(param.uuid, param.axis[0]["axis"])
            ax.text(x_min, legend_y - 0.35*i,
                    f"■ Axis {i+1}: {param.name}, "
                    f"{fmt_range(param.value + axis[0], param.value + axis[-1], param.unit)}"
                    f" ({len(axis)} pts)",
                    ha='left', va='top', fontsize=9, color=evo_colors[param.uuid])

        ax.set_xlim(x_min - 0.2, x_max)
        ax.set_ylim(y_min, baselines[0] + height + 0.6)
        ax.set_title(kwargs.get("title", sequence.name))
        ax.axis('off')

        return fig