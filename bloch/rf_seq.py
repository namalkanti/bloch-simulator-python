import numpy as np
import scipy as sp

def hard_pulses(b1, flip_angle, seperation, length, count, dt, gamma=26752):
    """
    Generates a series of hard rf pulses with a specific flip angle seperation and length.
    
    b1: Amplitude of rf pulse
    flip_angle: Array of flip angle for rf pulse in radians.
    seperation: Time seperation for different rf puslses in seconds.
    length: Total time for rf pulse length in seconds.
    count: Number of rf pulses in the sequence.
    dt: dt for each index
    gamma: Gamma value for environment default is set for hydrogen in rad/s/G.
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
