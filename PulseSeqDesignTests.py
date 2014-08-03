import unittest
import numpy as np
import scipy as sp

import scipy.io as sio

from bloch.pulse_seq_design import generate_readout_gradient, generate_phase_encode_gradient

from BlochTest import get_data_with_key

TEST_DIR = "test_data" 
TEST_FILE = "pulse"

g_max = 4
s_max = 15000
dt = 4e-6

class PulseSeqDesignTests(unittest.TestCase):

    def test_readout(self):
        """
        Compares readout generation to matlab results.
        """
        expected_gx = get_data_with_key(TEST_DIR, TEST_FILE, "gx")
        expected_rowin = get_data_with_key(TEST_DIR, TEST_FILE, "rowin")

        Nf = 64
        Fov_r  = 14
        bwpp = 1862.4
        gx, rowin = generate_readout_gradient(Nf, Fov_r, bwpp, g_max, s_max, dt)
        self.assertTrue(np.allclose(expected_gx, gx))
        self.assertTrue(np.allclose(expected_rowin, rowin))


    def test_phase(self):
        """
        Compares phase generation to matlab results.
        """
        expected_gpe = get_data_with_key(TEST_DIR, TEST_FILE, "gpe")
        expected_petable = get_data_with_key(TEST_DIR, TEST_FILE, "petable")

        Np = 32
        Fov_p = 7
        gpe, petable = generate_phase_encode_gradient(Np, Fov_p, g_max, s_max, dt)
        self.assertTrue(np.allclose(expected_gpe, gpe))
        self.assertTrue(np.allclose(expected_petable, petable))

if __name__ == "__main__":
    unittest.main()
