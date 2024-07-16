import matplotlib.pyplot as plt
import numpy as np
from math import pi, floor
from .units import *

class Shape:
    def __init__(
        self,
        shape: str,
        width_ms: float,
    ):
        # from sympy import Symbol, Function, diff, lambdify, sin
        import sympy as sp

        # Shape can be "Sine", "Damped_Sine", "Monophasic", "Square"
        # Only "Sine" and "Square" are supported currently
        self.shape = shape
        self.width_ms = width_ms
        if self.shape not in ["Sine", "Square"]:
            raise ValueError(
                "The only shapes supported currently are 'Sine' and 'Square'"
            )
        if self.width_ms <= 0:
            raise ValueError("Width must be > 0")
        t = sp.Symbol("t")
        f = sp.Function(t)
        if self.shape == "Sine":
            sin_freq_kHz = 2 * pi / width_ms # Angular frequency of sine wave
            f = -sp.sin(t * sin_freq_kHz) / sin_freq_kHz # Sine wave scaled so that derivative amplitude is 1
        elif self.shape == "Square":
            f = t # Derivative will be constant value of 1
        self.current_shape = sp.lambdify(t, f) # Define current function shape (agnostic to amplitude)
        self.efield_shape = sp.lambdify(t, sp.diff(f, t)) # Define electric field function shape (amplitude of 1; scaling comes later)

from typing import Union

class Pattern:
    def __init__(
        self,
        pshape: Shape,  # Shape of the pulse
        pattern: Union[str, None] = None,  # Pattern of the set
        npulses: Union[int, None] = None,  # Number of pulses in a set
        interval_ms: Union[float, None] = None,  # Duration of interval between pulses in a set
        p_onset_interval_ms: Union[float, None] = None,  # Duration of interval between onset of pulses in a set
        set_freq_Hz: Union[float, None] = None,  # Frequency of pulse onsets in a set
    ):
        # ... rest of your code ...

# class Pattern:
#     def __init__(
#         self,
#         pshape: Shape,  # Shape of the pulse
#         pattern: str | None = None,  # Pattern of the set
#         npulses: int | None = None,  # Number of pulses in a set
#         interval_ms: float | None = None,  # Duration of interval between pulses in a set
#         p_onset_interval_ms: float | None = None,  # Duration of interval between onset of pulses in a set
#         set_freq_Hz: float | None = None,  # Frequency of pulse onsets in a set
#     ):
        """If pattern is "Single" then npulses set to 1 (then interval_ms, p_onset_interval_ms, & set_freq_Hz are meaningless)
        If "TBS" then npulses set to 3 (Theta burst stimulation)

        interval_ms = p_onset_interval_ms - pshape.width_ms
        p_onset_interval_ms = 1/set_freq_Hz
        Highest priority when defined | interval_ms > p_onset_interval_ms > set_freq_Hz | lowest priority"""

        # Set attributes
        self.pshape = pshape
        self.pattern = pattern
        width_ms = self.pshape.width_ms

        # Implement pattern
        if self.pattern == "Single":
            interval_ms = 0 # Meaningless for this pattern, but marks it as set
            npulses = 1
        elif self.pattern == "TBS":
            # An interval parameter must still be defined
            npulses = 3
        
        # Check that npulses and pattern are valid
        if type(npulses) != int or npulses < 1:
            raise ValueError(
                "npulses must be defined as a positive non-zero integer or pulse pattern must be categorized as either 'Single' or 'TBS'"
                )
        self.npulses = npulses

        # Check that interval_ms, p_onset_interval_ms, and set_freq_Hz are valid
        if interval_ms == None:
            if p_onset_interval_ms == None:
                if set_freq_Hz == None:
                    raise ValueError(
                        "Pulse interval, onset interval, or set frequency of pattern must be defined"
                    )
                else:
                    if set_freq_Hz > 1/width_ms * kHz: # Comparison in Hz
                        raise ValueError(
                            f"Set freq [{set_freq_Hz} Hz] must be <= 1/(pulse width) [{1/width_ms * kHz} Hz]"
                        )
                    interval_ms = 1/set_freq_Hz * s - width_ms
            else:
                if p_onset_interval_ms < width_ms:
                    raise ValueError(
                        f"Pulse onset interval [{p_onset_interval_ms} ms] must be >= pulse width [{width_ms} ms]"
                    )
                interval_ms = p_onset_interval_ms - width_ms
        if interval_ms < 0:
            raise ValueError(f"Pulse interval of pattern [{interval_ms} ms] must be >= 0 ms")
        
        # Set remaining attributes
        self.interval_ms = interval_ms
        self.p_onset_interval_ms = self.interval_ms + width_ms
        self.set_freq_Hz = 1 / self.p_onset_interval_ms * kHz


def efield(
    freq_Hz: float,  # Frequency of pulse sets in Hz
    duration_ms: float,  # Duration of waveform
    pulse_resolution_ms: float,  # Duration of time step
    stim_start_ms: float,  # Initial waiting period
    stim_end_ms: float, # Time when stimulation ends
    ef_amp_mV_per_um: float,  # Amplitude of the max E-field in the desired waveform
    pat: Pattern,  # Pattern object containing data on the waveform
):
    rd = 9 # Rounding precision to correct floating point error (rounds to 10^-rd ms)
    pulsef = pat.pshape.efield_shape  # E-field pulse function (SymPy)
    pwidth_ms = pat.pshape.width_ms  # Duration of one pulse
    p_onset_interval_ms = pat.p_onset_interval_ms # Duration of interval between the onset of pulses in a set
    inter_p_interval_ms = pat.interval_ms # Duration of interval between pulses within a set
    npulses = pat.npulses # Number of pulses in one set
    set_width_ms = (npulses - 1) * inter_p_interval_ms + npulses * pwidth_ms # Duration of one set of pulses (time of intervals + time of pulses)
    set_onset_interval_ms = 1 / freq_Hz * s # Duration of interval between the onset of sets of pulses
    inter_set_interval_ms = set_onset_interval_ms - set_width_ms # Duration of interval between sets of pulses

    if stim_end_ms == None:
        stim_end_ms = duration_ms

    # Check that inputs are valid
    if freq_Hz <= 0:
        raise ValueError("Frequency must be > 0")
    if duration_ms <= 0:
        raise ValueError("Duration must be > 0")
    if pulse_resolution_ms <= 0:
        raise ValueError("pulse_resolution must be > 0")
    if inter_set_interval_ms < 0:
        raise ValueError(
            f"Duration of pulse set [{set_width_ms} ms] must be <= interval between pulse onset (1/frequency) [{set_onset_interval_ms} ms]"
        )

    # Construct waveform of electric field, taking advantage of linear interpolation between time points

    # Initialize variables for waveform construction
    time = [] # Points in time
    wav = []  # E-field waveform at each point in time

    if stim_start_ms > stim_end_ms:
        # Waveform is silent for entire duration
        time = [0, duration_ms]
        wav = [0, 0]
        return [wav, time]

    if stim_start_ms > 0:
        # Start of silent period
        time.append(0)
        wav.append(0)

    nstep_p = int(round(pwidth_ms / pulse_resolution_ms))  # Number of time steps within a pulse
    npoint_p = nstep_p + 1  # Number of time points within a pulse
    pulse_t = np.linspace(0, pwidth_ms, npoint_p) # Time points within a pulse starting at t=0
    pulse = np.array([pulsef(t) for t in pulse_t]) * ef_amp_mV_per_um # Sample points of the pulse waveform; scale by ef_amp
    
    cur_t = stim_start_ms
    pulse_start_times = [] # List of pulse start times
    while cur_t < stim_end_ms: # Iterate through duration and build pulse_start_times
        for pcount in range(npulses):
            pulse_start_times.append(cur_t)
            if pcount == npulses-1: # If at the end of the set
                cur_t += inter_set_interval_ms + pwidth_ms
            else: # If in the middle of a set
                cur_t += p_onset_interval_ms

    pulse_end_time = 0
    for i, pulse_start in enumerate(pulse_start_times):
        # Pre-pulse buffer
        last_silent_t = pulse_start - pulse_resolution_ms # Last time point of silent period before pulse start
        if last_silent_t > pulse_end_time: # If we've advanced far enough to need to specify the end of the silent period
            # Write end of silent period
            time.append(last_silent_t)
            wav.append(0)

        # Write pulse
        time.extend(pulse_t + pulse_start)
        wav.extend(pulse)

        # Post-pulse buffer
        pulse_end_time = pulse_start + pwidth_ms # Time of the end of the pulse
        start_silent_t = pulse_end_time + pulse_resolution_ms # First time point of silent period after pulse end
        if (i == len(pulse_start_times)-1 # If this is the last pulse
            or start_silent_t < pulse_start_times[i+1]): # or if we will advance far enough to need to specify the start of a silent period
            # Write start of silent period
            time.append(start_silent_t)
            wav.append(0)
    # Buffers necessary for simulation time steps when outside of a pulse

    if time[-1] < duration_ms: # If the time course does not last the full duration
        # Place a silent period until the end of the simulation
        time.append(duration_ms)
        wav.append(0)
    time = [round(t, rd) for t in time] # Correct floating point error

    return [wav, time]


def get_efield(
    freq_Hz: float,
    duration_ms: float,
    pulse_resolution_ms: float = 1e-3,
    stim_start_ms: float = 0.,
    stim_end_ms: Union[float, None] = None,
    ef_amp_V_per_m: float = 100.,
    pshape: str = "Sine",
    width_ms: float = 100e-3,
    pat: str = "Single",
    npulses: int = None,
    interval_ms: float = None,
    p_onset_interval_ms: float = None,
    set_freq_Hz: float = None,
    plot: bool = False,
):
    """
    freq_Hz: Frequency of TMS pulses in Hz
    duration_ms: Duration of simulation in ms
    pulse_resolution_ms: Temporal resolution of pulses in ms (independent of simulation dt)
    stim_start_ms: Time of first pulse in ms
    stim_end_ms: Time when stimulation ends in ms
    ef_amp_V_per_m: Amplitude of pulse in V/m
    pshape: Qualitative description of waveform (see Shape class)
    width_ms: Period of waveform in ms
    pat: Qualitative description of stimulation pattern (see Pattern class)
    npulses: Number of pulses in one set of a pattern
    interval_ms: Duration of interval between pulses in a set in ms
    p_onset_interval_ms: Duration of interval between onset of pulses in a set in ms
    set_freq_Hz: Frequency of pulse onsets in a set in Hz
    plot: Whether to plot the results
    
    Returns waveform in mV/um (or V/mm)

    Returns time course in ms
    """

    wav, time = efield(
        freq_Hz=freq_Hz,
        duration_ms=duration_ms,  
        pulse_resolution_ms=pulse_resolution_ms, 
        stim_start_ms=stim_start_ms,  
        stim_end_ms=stim_end_ms,
        ef_amp_mV_per_um=ef_amp_V_per_m / (mV/um), # Convert from V/m to mV/um (or V/mm)
        pat=Pattern(
            pshape=Shape(shape=pshape, width_ms=width_ms),
            pattern=pat,
            npulses=npulses,
            interval_ms=interval_ms,
            p_onset_interval_ms=p_onset_interval_ms,
            set_freq_Hz=set_freq_Hz,
        ),
    )

    if plot:
        plt.figure()
        plt.step(time, wav, where="post")
        plt.xlabel("Time (ms)")
        plt.ylabel("E-field (mV/um)")
        plt.show()

    return wav, time


# Testing
def main():
    get_efield(
        freq_Hz=10,
        duration_ms=500,
        width_ms=100e-3,
        pshape="Sine",
        pat="Single",
        # npulses=2,
        # interval_ms=100e-3,
    )

    plt.show()

    """,
    dt_ms=25e-3,
    pat="Unique",
    npulses=2,
    interval=100e-6,
    p_onset_interval=None,
    set_freq_Hz=None"""

    """plt.step(time, wav, where='post')
    plt.show()"""


if __name__ == "__main__":
    main()
