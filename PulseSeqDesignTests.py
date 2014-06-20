import unittest
import numpy as np
import scipy as sp

import scipy.io as sio

from bloch.pulse_seq_design import genReadoutGradient, genPEGradient

g_max = 4
s_max = 15000
dt = 4e-6

class PulseSeqDesignTests(unittest.TestCase):

    def test_readout(self):
        """
        Compares readout generation to matlab results.
        """
        expected_gx = sio.loadmat("bloch/test_data/pulse/gx.mat")["gro"].ravel()
        expected_rowin = sio.loadmat("bloch/test_data/pulse/rowin.mat")["rowin"].ravel()

        Nf = 64
        Fov_r  = 14
        bwpp = 1862.4
        gx, rowin = genReadoutGradient(Nf, Fov_r, bwpp, g_max, s_max, dt)
        self.assertTrue(np.allclose(expected_gx, gx))
        self.assertTrue(np.allclose(expected_rowin, rowin))


    def test_phase(self):
        """
        Compares phase generation to matlab results.
        """
        expected_gpe = sio.loadmat("bloch/test_data/pulse/gpe.mat")["gpe"].ravel()
        expected_petable = sio.loadmat("bloch/test_data/pulse/petable.mat")["petable"].ravel()

        Np = 32
        Fov_p = 7
        gpe, petable = genPEGradient(Np, Fov_p, g_max, s_max, dt)
        self.assertTrue(np.allclose(expected_gpe, gpe))
        self.assertTrue(np.allclose(expected_petable, petable))

if __name__ == "__main__":
    unittest.main()
