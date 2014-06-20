import numpy as np
import scipy as sp

from bloch.minTimeGradientArea import minTimeGradientArea

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
    
    def generate_readout(self, Nf, Fov_r, bwpp):
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
        gpre = minTimeGradientArea(areaTrapz/2, self._g_max, self._s_max, self._dt)

        rowin = gpre.size + 1 + np.asarray(idx2)

        gro = np.concatenate((-gpre, np.array([0]), gro))

        return gro, rowin

    def generate_phase_encodes(self, Np, For_p):
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
        grpe = minTimeGradientArea(area, self._g_max, self._s_max, self._dt)

        petable = np.arange(Np/2 - .5, -Np/2 + .5 - 1, -1) / (Np/2) 

        return grpe, petable

def genReadoutGradient(Nf, FOVr, bwpp, Gmax, Smax, dt):
    """
    Input:
    Nf: The number of frequency encodes
    FOVr: The desired field of view (cm)
    bwpp: The desired bandwidth per pixel (Hz/pixel)
    Gmax: The maximum gradient (G/cm)
    Smax: The maximum slew rate (G/cm/s)
    dt: The duration of each sample (s)

    Output:
    gro: Gradient waveform.
    rowin: Indices corresponding to the readout portion of the gradient.
    """
    gradient_generator = GradientCalculator(Gmax, Smax, dt)
    return gradient_generator.generate_readout(Nf, FOVr, bwpp)

def genPEGradient(Np, FOVp, Gmax, Smax, dt):
    """
    Input:
    Np: The number of phase encodes
    FOVp: The desired field of view (cm)
    Gmax: The maximum gradient (G/cm)
    Smax: The maximum slew rate (G/cm/s)
    dt: The duration of each sample (s)

    Output:
    grpe: Gradient waveform.
    petable: Indices corresponding to the readout portion of the gradient.
    """
    gradient_generator = GradientCalculator(Gmax, Smax, dt)
    return gradient_generator.generate_phase_encodes(Np, FOVp)
