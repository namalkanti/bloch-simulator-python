import unittest 

import numpy as np
import scipy.io as sio

from bloch import bloch
from pulse_seq_design import genReadoutGradient
from pulse_seq_design import genPEGradient

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

    def test_hw5_bloch_sim_first_part(self):
        matlab = sio.loadmat(TEST_DIR.format("hw5", "hw5_img"))

        dp = matlab["dp"]
        mx_0 = matlab["mx"].ravel()
        my_0 = matlab["my"].ravel()
        mz_0 = matlab["mz"].ravel()

        expected_mx1 = sio.loadmat(TEST_DIR.format("hw5", "mxa.mat"))["mx1"].ravel()
        expected_my1 = sio.loadmat(TEST_DIR.format("hw5", "mya.mat"))["my1"].ravel()
        expected_mz1 = sio.loadmat(TEST_DIR.format("hw5", "mza.mat"))["mz1"].ravel()

        expected_mx2 = np.transpose(sio.loadmat(TEST_DIR.format("hw5", "mxb.mat"))["mx1"])
        expected_my2= np.transpose(sio.loadmat(TEST_DIR.format("hw5", "myb.mat"))["my1"])
        expected_mz2= np.transpose(sio.loadmat(TEST_DIR.format("hw5", "mzb.mat"))["mz1"])

        Nf = 64
        Np = 32
        Nrf = 92
        Fov_r = 14
        Fov_p = 7
        g_max = 4
        s_max = 15000
        dt = 4e-6
        bwpp = 1862.4
        gamma = 4257
        flip = 90

        gx, rowin = genReadoutGradient(Nf, Fov_r, bwpp, g_max, s_max, dt)
        gpe, petable = genPEGradient(Np, Fov_p, g_max, s_max, dt)
        rf_90 = np.ones(Nrf) * (flip/360) / (Nrf * dt * gamma)
        gy = np.zeros(gx.size)
        gy[:gpe.size] = gpe

        g = np.asarray([gx, gy * -petable[0]])
        mx1, my1, mz1 = bloch(rf_90, rf_90 * 0, dt, 100, 100, 0, dp, 0, mx_0, my_0, mz_0)
        mx2, my2, mz2 = bloch(gx * 0, g, dt, 100, 100, 0, dp, 2, mx1, my1, mz1)

        self.assertTrue(np.allclose(expected_mx1, mx1))
        self.assertTrue(np.allclose(expected_my1, my1))
        self.assertTrue(np.allclose(expected_mz1, mz1))

        self.assertTrue(np.allclose(expected_mx2, mx2))
        self.assertTrue(np.allclose(expected_my2, my2))
        self.assertTrue(np.allclose(expected_mz2, mz2))

if __name__ == "__main__":
    unittest.main()
