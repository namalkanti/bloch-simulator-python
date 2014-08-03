import unittest

import numpy as np
import scipy as sp

from bloch.rf_seq import hard_pulses

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

if __name__ == "__main__":
    unittest.main()
