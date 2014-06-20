import unittest

import numpy as np
import scipy as sp

import scipy.io as sio

from minTimeGradientArea import minTimeGradientArea

class GradientAreaTest(unittest.TestCase):

    def test_waveforms(self):
        """
        Tests waveforms against matlab code results.
        """
        expected_gy = sio.loadmat("test_data/gradient/gy_sample.mat")["gy"].ravel()

        Np = 32
        fov = 7
        smax = 15000
        gmax = 4
        gamma = 4257
        dt = 4e-6

        kmax = 1/(fov/Np)/2
        area = kmax / gamma
        gy = minTimeGradientArea(area, gmax, smax, dt)
        self.assertTrue(np.allclose(expected_gy, gy))

if __name__ == "__main__":
    unittest.main()
