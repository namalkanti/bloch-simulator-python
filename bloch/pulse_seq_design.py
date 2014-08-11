import numpy as np
import scipy as sp

from bloch.min_time_gradient import minimum_time_gradient

class GradientCalculator():
    """
    Class used to generate readout and phase gradients for an MRI simulation.
    """

    def __init__(self, g_max, s_max, dt):
        """
        Input:
        g_max: The maximum gradient (G/cm)
        s_max: The maximum slew rate (G/cm/s)
        dt: The duration of each sample (s)
        """
        self._g_max = g_max
        self._s_max = s_max
        self._dt = dt
        self._gamma = 4257
    
    def generate_readout(self, Nf, Fov_r, bwpp, negate=-1, pre_delay=0, readout_delay=1):
        """
        Input:
        Nf: The number of frequency encodes
        Fov_r: The desired field of view (cm)
        bwpp: The desired bandwidth per pixel (Hz/pixel)
        
        Output:
        gro: Gradient waveform.
        rowin: Indices corresponding to the readout portion of the gradient.
        """
        res = Fov_r / Nf
        Wkx = 1 / res
        area = Wkx / self._gamma

        G = bwpp / res / self._gamma
        Tro = Wkx / self._gamma / G
        Tramp = G / self._s_max

        t1 = Tramp
        t2 = t1 + Tro
        T = Tramp * 2 + Tro

        N = int(np.floor(T / self._dt))
        t = np.arange(1, N + 1) * self._dt

        idx1 = np.where(t <  t1) 
        idx2 = np.where((t >= t1) & (t < t2))
        idx3 = np.where(t >= t2)

        gro = np.zeros(N)
        gro[idx1] = self._s_max * t[idx1]
        gro[idx2] = G
        gro[idx3] = T * self._s_max - self._s_max * t[idx3]

        areaTrapz = (T + Tro) * G/2
        gpre = minimum_time_gradient(areaTrapz/2, self._g_max, self._s_max, self._dt)

        rowin = pre_delay + gpre.size + readout_delay + np.asarray(idx2)

        gro = np.concatenate((np.zeros(pre_delay), negate * gpre, np.zeros(readout_delay), gro))

        return gro, rowin

    def generate_spin_echo_readout(self, Nf, Fov_r, bwpp, rf_duration, te):
        """
        Generates a spin echo readout gradient with specified te as time echo.
        RF duration and te should be provided in seconds
        Input:
        Nf: The number of frequency encodes
        Fov_r: The desired field of view (cm)
        bwpp: The desired bandwidth per pixel (Hz/pixel)
        rf_duration: Duration of rf pulse(delays prewinder)
        te: Echo time for readout gradient.
        
        Output:
        gro: Gradient waveform.
        rowin: Indices corresponding to the readout portion of the gradient.
        """
        pre_delay = self._dt * rf_duration 
        readout_delay = self._dt * te
        return self.generate_readout(Nf, Fov_r, bwpp, 1, pre_delay, readout_delay)

    def generate_phase_encodes(self, Np, For_p, delay=0):
        """
        Input:
        Np: The number of phase encodes
        Fov_p: The desired field of view (cm)

        Output:
        grpe: Gradient waveform.
        petable: Indices corresponding to the readout portion of the gradient.
        """

        kmax = 1 / (For_p / Np) / 2
        area = kmax / self._gamma
        grpe = minimum_time_gradient(area, self._g_max, self._s_max, self._dt)

        petable = delay + np.arange(Np/2 - .5, -Np/2 + .5 - 1, -1) / (Np/2) 

        grpe = np.concatenate((np.zeros(delay), grpe))


        return grpe, petable

    def generate_delayed_phase_encodes(self, Np, For_p, delay_time):
        """
        Accepts a delay time in seconds. Otherwise the same as generate_phase_encodes.
        Input:
        Np: The number of phase encodes
        Fov_p: The desired field of view (cm)
        delay_time: Time is seconds until phase encodes start.

        Output:
        grpe: Gradient waveform.
        petable: Indices corresponding to the readout portion of the gradient.
        """
        delay_samples = self._dt * delay_time
        return self.generate_phase_encodes(Np, For_p, delay_samples)


def generate_readout_gradient(Nf, fov_r, bwpp, g_max, s_max, dt):
    """
    Input:
    Nf: The number of frequency encodes
    fov_r: The desired field of view (cm)
    bwpp: The desired bandwidth per pixel (Hz/pixel)
    g_max: The maximum gradient (G/cm)
    s_max: The maximum slew rate (G/cm/s)
    dt: The duration of each sample (s)

    Output:
    gro: Gradient waveform.
    rowin: Indices corresponding to the readout portion of the gradient.
    """
    gradient_generator = GradientCalculator(g_max, s_max, dt)
    return gradient_generator.generate_readout(Nf, fov_r, bwpp)

def generate_phase_encode_gradient(Np, fov_p, g_max, s_max, dt):
    """
    Input:
    Np: The number of phase encodes
    fov_p: The desired field of view (cm)
    g_max: The maximum gradient (G/cm)
    s_max: The maximum slew rate (G/cm/s)
    dt: The duration of each sample (s)

    Output:
    grpe: Gradient waveform.
    petable: Indices corresponding to the readout portion of the gradient.
    """
    gradient_generator = GradientCalculator(g_max, s_max, dt)
    return gradient_generator.generate_phase_encodes(Np, fov_p)
