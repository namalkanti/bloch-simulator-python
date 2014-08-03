import os.path
import unittest

import numpy as np
import scipy as sp

import scipy.io as sio

from bloch.min_time_gradient import minimum_time_gradient

from BlochTest import get_data_with_key

TEST_DIR = "test_data"
TEST_FILE = "gradient"

class GradientAreaTest(unittest.TestCase):

    def test_waveforms(self):
        """
        Tests waveforms against matlab code results.
        """
        expected_gy = get_data_with_key(TEST_DIR, TEST_FILE, "gy") 

        Np = 32
        fov = 7
        smax = 15000
        gmax = 4
        gamma = 4257
        dt = 4e-6

        kmax = 1/(fov/Np)/2
        area = kmax / gamma
        gy = minimum_time_gradient(area, gmax, smax, dt)
        self.assertTrue(np.allclose(expected_gy, gy))

if __name__ == "__main__":
    unittest.main()
