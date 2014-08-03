import os.path
import unittest 

import numpy as np
import scipy.io as sio

from bloch.bloch import bloch

TEST_DIR = "test_data"

class BlochTest(unittest.TestCase):
    """
    Tests basic functionality of Python
    Bloch simulator. All test results are based 
    on results of Matlab function. Any errors in that
    evaluation will be replicated here. 
    """

    def test_bloch_sim_demo(self):
        """
        Runs a simple Bloch simulator run.
        """
        file_name = "basic_bloch"
        mx_demo = get_data_with_key(TEST_DIR, file_name, "mx_demo")
        my_demo = get_data_with_key(TEST_DIR, file_name, "my_demo")
        mz_demo = get_data_with_key(TEST_DIR, file_name, "mz_demo")

        b1 = get_data_with_key(TEST_DIR, file_name, "b1_demo")
        g = get_data_with_key(TEST_DIR, file_name, "g_demo")

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

    def test_bloch_sim_a(self):
        """
        Runs another simple Bloch simulator run.
        """
        mx_a = get_data_with_key(TEST_DIR, "basic_bloch", "mx_a")
        my_b = get_data_with_key(TEST_DIR, "basic_bloch", "my_a")
        mz_c = get_data_with_key(TEST_DIR, "basic_bloch", "mz_a")

        b1 = get_data_with_key(TEST_DIR, "basic_bloch", "b1_a")
        g = get_data_with_key(TEST_DIR, "basic_bloch", "g_a")

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
    
    def test_ssfptransiest(self):
        """
        Runs an SSFP response calculation using bloch.m
        """
        expected_mxss = get_data_with_key(TEST_DIR, "ssfptransient", "mxss")
        expected_myss = get_data_with_key(TEST_DIR, "ssfptransient", "myss")
        expected_mzss = get_data_with_key(TEST_DIR, "ssfptransient", "mzss")

        expected_mx = get_data_with_key(TEST_DIR, "ssfptransient", "mx")
        expected_my = get_data_with_key(TEST_DIR, "ssfptransient", "my")
        expected_mz = get_data_with_key(TEST_DIR, "ssfptransient", "mz")

        TR = .005 #Seconds
        Trf = 0.0001 #100 us "hard" RF pulse
        alpha = 60 #Degrees
        gamma = 4258 #Hz/G
        T1 = 1 #Seconds
        T2 = .2 #Seconds
        freq = np.arange(-200, 201) #Hz
        N = 100
        Tpad = (TR - Trf)/2 #Seconds

        t = np.asarray([Tpad, Trf, Tpad])
        b1 = np.concatenate((np.zeros(1), np.asarray([np.pi/180*alpha/Trf/gamma/2/np.pi]), np.zeros(1)))

        mxss, myss, mzss = bloch(b1, 0*b1, t, T1, T2, freq, 0, 1)

        mx, my, mz = bloch(1.0j * np.max(b1)/2, 0, Trf, T1, T2, freq, 0, 0)

        for _ in range(N):
            mx, my, mz = bloch(b1, 0*b1, t, T1, T2, freq, 0, 0, mx, my, mz)

        self.assertTrue(np.allclose(expected_mxss, mxss))
        self.assertTrue(np.allclose(expected_myss, myss))
        self.assertTrue(np.allclose(expected_mzss, mzss))

        self.assertTrue(np.allclose(expected_mx, mx))
        self.assertTrue(np.allclose(expected_my, my))
        self.assertTrue(np.allclose(expected_mz, mz))

def get_data(directory, file_name):
    """
    Grabs data structure from matlab data file.
    """
    file_name = "{0}.npz".format(file_name)
    data = np.load(os.path.join(directory, file_name))
    return data

def get_data_with_key(directory, file_name, key):
    """
    Grabs data structure from matlab data file and key.
    """
    data = get_data(directory, file_name)[key]
    return np.transpose(data)

if __name__ == "__main__":
    unittest.main()
