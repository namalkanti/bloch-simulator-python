import numpy as np
import scipy as sp

import matplotlib.pyplot as plt
import scipy.io as sio

from numpy.random import randn

from bloch import bloch

TEST_DIR = "test_data/{0}/{1}"

def main():
    mx_demo = sio.loadmat(TEST_DIR.format("hw4", "mx_demo"))["mx"].ravel()
    my_demo = sio.loadmat(TEST_DIR.format("hw4", "my_demo"))["my"].ravel()
    mz_demo = sio.loadmat(TEST_DIR.format("hw4", "mz_demo"))["mz"].ravel()

    b1 = sio.loadmat(TEST_DIR.format("hw4", "b1_demo"))["b1"]
    g = sio.loadmat(TEST_DIR.format("hw4", "g_demo"))["g"]

    dt = 4e-6
    t1 = 100e-3
    t2 = 50e-3
    df = 0
    dp = 0
    mode = 2

    mx_0 = 0
    my_0 = 0
    mz_0 = 1

    mx, my, mz = bloch(b1, g, dt, t1, t2, df, dp, mode, mx_0, my_0, mz_0)
    x_values = np.arange(0, 100)
    mz = np.exp(x_values)
    plt.plot(x_values, abs(mz))
    plt.show()

if __name__ == "__main__":
    main()
