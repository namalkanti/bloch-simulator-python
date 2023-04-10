import numpy as np
import scipy as sp
import scipy.signal as sig

def hard_pulses(b1, flip_angle, seperation, length, count, dt, gamma=26752):
    """
    Generates a series of hard rf pulses with a specific flip angle seperation and length.
    
    Input:
    b1: Amplitude of rf pulse
    flip_angle: Array of flip angle for rf pulse in radians.
    seperation: Time seperation for different rf puslses in seconds.
    length: Total time for rf pulse length in seconds. If length is too short, an error will be thrown.
    count: Number of rf pulses in the sequence.
    dt: dt for each index
    gamma: Gamma value for environment default is set for hydrogen in rad/s/G.

    Output:
    rf: Hard pulse rf signal
    """
    total_units = int(length / dt)
    rf = np.zeros(total_units)
    seperation_units = int(seperation / dt)
    for val in range(count):
        pulse_time = flip_angle[val] / (gamma * b1)
        pulse_units = int(pulse_time / dt)
        start_idx = val * (seperation_units)
        rf[start_idx:start_idx+pulse_units] = b1 
    return rf

def sinc_pulse(timebandwidth, flip_angle, duration, dt, gamma=26747.52):
    """
    Generates an rf pulse with specified tbw, duration, and dt.

    Input:
    timebandwidth: Timebandwidth of desired rf pulse
    flip_angle: Flip angle of desired pulse
    duration: RF pulse duration in seconds.
    dt: dt value for rf samples
    gamma: Gamma value. Defaults to 2 * pi * 4257(hydrogren default in radians)

    Output:
    rf: Sinc pulse rf signal
    """
    samples = int(duration / dt)
    theta = np.linspace(-timebandwidth/2, timebandwidth/2, samples+2)
    rf = np.sinc(theta[1:-1]) * sig.hann(samples)
    rf = flip_angle * (rf/np.sum(rf))
    rf /= (gamma * dt)
    return rf
