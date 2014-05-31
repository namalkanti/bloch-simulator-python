import unittest 

import numpy as np
import scipy.io as sio

from bloch import bloch

TEST_DIR = "test_data/{0}/{1}"

class BlochTest(unittest.TestCase):
    """
    Tests basic functionality of Python
    Bloch simulator. All test results are based 
    on results of Matlab function. Any errors in that
    evaluation will be replicated here. 
    """

    def test_hw4_bloch_sim_demo(self):
        mx_demo = sio.loadmat(TEST_DIR.format("hw4", "mx_demo"))["mx"].ravel()
        my_demo = sio.loadmat(TEST_DIR.format("hw4", "my_demo"))["my"].ravel()
        mz_demo = sio.loadmat(TEST_DIR.format("hw4", "mz_demo"))["mz"].ravel()

        b1 = sio.loadmat(TEST_DIR.format("hw4", "b1_demo"))["b1"].ravel()
        g = np.transpose(sio.loadmat(TEST_DIR.format("hw4", "g_demo"))["g"]).ravel()

        dt = 4e-6
        t1 = 100e-3
        t2 = 50e-3
        df = 0
        dp = 0
        mode = 2

        mx_0 = 0
        my_0 = 0
        mz_0 = 1

        mx, my, mz = bloch(b1, g, dt, t1, t2, df, dp, mode, mx_0, mx_0, mx_0)
        self.assertTrue(np.allclose(mx_demo, mx));
        self.assertTrue(np.allclose(my_demo, my));
        self.assertTrue(np.allclose(mz_demo, mz));

    def test_hw4_bloch_sim_a(self):
        mx_a = sio.loadmat(TEST_DIR.format("hw4", "mx_a"))["mx"].ravel()
        my_b = sio.loadmat(TEST_DIR.format("hw4", "my_a"))["my"].ravel()
        mz_c = sio.loadmat(TEST_DIR.format("hw4", "mz_a"))["mz"].ravel()

        b1 = sio.loadmat(TEST_DIR.format("hw4", "B1_a"))["B1"].ravel()
        g = np.transpose(sio.loadmat(TEST_DIR.format("hw4", "G_a"))["G"]).ravel()

        dt = 4e-6
        t1 = 30e-3
        t2 = 15e-3
        df = 0
        dp = .2
        mode = 2

        mx_0 = 0
        my_0 = 0
        mz_0 = 1

        mx, my, mz = bloch(b1, g, dt, t1, t2, df, dp, mode, mx_0, my_0, mz_0)
        self.assertTrue(np.allclose(mx_a, mx));
        self.assertTrue(np.allclose(my_b, my));
        self.assertTrue(np.allclose(mz_c, mz));

if __name__ == "__main__":
    unittest.main()
