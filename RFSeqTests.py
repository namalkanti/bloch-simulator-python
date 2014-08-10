import unittest

import numpy as np
import scipy as sp

from bloch.bloch import bloch
from bloch.rf_seq import hard_pulses, sinc_pulse

from BlochTest import get_data_with_key

TEST_DIR = "test_data"

class RFSeqTests(unittest.TestCase):
    """
    Tests the rf sequences included with the module.
    """

    def test_hard_pulses(self):
        """
        Tests generation of a series of hard pulses.
        """

        gamma = 26752
        dt = 4e-6
        b1_max = .16
        expected_pulses = np.zeros(int(100e-3 / dt))
        pulse_time = (np.pi / 2) / (gamma * .16)  
        pulse_length = int(pulse_time / dt)
        expected_pulses[:pulse_length] = b1_max
        second_pulse = int(20e-3/dt)
        expected_pulses[second_pulse:second_pulse+pulse_length] = b1_max

        b1 = .16
        flip_angle = np.ones(2) * np.pi / 2
        spin_echo_pulses = hard_pulses(b1, flip_angle, 20e-3, 100e-3, 2, 4e-6)
        self.assertTrue(np.allclose(expected_pulses, spin_echo_pulses))

    def test_sinc_pulse(self):
        """
        Tests generation of sinc pulse.
        """
        expected_mxy = get_data_with_key(TEST_DIR, "rf_test", "mxy")

        duration = 3.2e-3
        dt = duration / 256
        flip = np.pi / 2
        tbw = 8

        rf = sinc_pulse(tbw, flip, duration, dt)
        g = np.ones(rf.size) * .6
        dp = np.linspace(-2, 2, 512)
        mx0 = np.zeros(512)
        my0 = np.zeros(512)
        mz0 = np.ones(512)

        mx, my, mz = bloch(rf, g, dt, 100, 100, 0, dp, 0, mx0, my0, mz0)
        mxy = mx + 1.0j*my
        np.allclose(expected_mxy, mxy)



if __name__ == "__main__":
    unittest.main()
