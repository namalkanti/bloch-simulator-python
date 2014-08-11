import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

class GradientMinimumTimeEstimator():
    """
    Class to generate gradient waveforms.
    """

    def __init__(self, area, G, S, dt):
        """
        Input:
        Area: Desired area for gradient waveform.
        G: Maximum strenght of gradient.
        S: Maximum slew rate for gradient.
        dt:Sampling interval
        """
        self._area = float(area)
        self._max_gradient = float(G)
        self._slew_rate = float(S)
        self._dt = dt
        self.calculate_waveform()

    def calculate_waveform(self):
        """
        Calculates waveform and assignes to instance variable.
        """
        triangle_area = self._max_gradient**2 / self._slew_rate

        if self._area <= triangle_area:
            t1 = np.sqrt(self._area / self._slew_rate)
            T = 2 * t1
            N = int(np.floor(T/self._dt))
            t = np.arange(1, N+1) * self._dt

            idx1 = np.where(t < t1)
            idx2 = np.where(t >= t1)

            g = np.zeros(N)
            g[idx1] = self._slew_rate * t[idx1]
            g[idx2] = 2 * np.sqrt(self._area * self._slew_rate) - self._slew_rate * t[idx2]

        else:
            t1 = self._max_gradient / self._slew_rate
            t2 = self._area / self._max_gradient
            t3 = self._area / self._max_gradient + (self._max_gradient / self._slew_rate)

            T = t3
            N = int(np.floor(T/self._dt))
            t = np.arange(1, N+1) * self._dt

            idx1 = np.where(t < t1)
            idx2 = np.where((t>=t1) & (t < t2))
            idx3 = np.where(t>=t2)

            g = np.zeros(N)
            g[idx1] = self._slew_rate * t[idx1]
            g[idx2] = self._max_gradient
            g[idx3] = (self._area/self._max_gradient + self._max_gradient/self._slew_rate) * self._slew_rate - self._slew_rate * t[idx3]

        self._waveform = g

    def plot(self):
        """
        Plots gradient waveform.
        """
        plt.plot(self._time, self._waveform)

    def get_waveform(self):
        """
        Returns gradient waveform.
        """
        return self._waveform

def minimum_time_for_area(area, Gmax, Smax, dt):
    """
    Returns the minumum time to achieve the desired gradient area with the given conditions.
    """
    estimator = GradientMinimumTimeEstimator(area, Gmax, Smax, dt)
    return estimator.get_total_time()

def minimum_time_gradient(area, Gmax, Smax, dt):
    """
    Returns minimum time gradient waveform for given parameters.
    """
    estimator = GradientMinimumTimeEstimator(area, Gmax, Smax, dt)
    return estimator.get_waveform()

def main():
    g1 = GradientMinimumTimeEstimator(6e-4, 4, 15000, 4e-6)
    g2 = GradientMinimumTimeEstimator(6e-4, 1, 5000, 4e-6)
    plt.figure()
    g1.plot()
    g2.plot()

if __name__ == "__main__":
    main()
    plt.show()
