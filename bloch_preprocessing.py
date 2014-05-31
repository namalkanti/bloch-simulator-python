import numpy as np
import scipy as sp

#Functions to handle preprocessing for bloch simulator arguments.

def process_gradient_argument(gr):
    """
    Takes in a gradient argument and returns directional gradients.
    If gradients don't exist, returns array of zeros.
    """
    if 1 == len(gr.shape):
        return gr, np.zeros(gr.size), np.zeros(gr.size)
    elif 3 == gr.shape[0]:
        return gr[0], gr[1], gr[2]
    elif 2 == gr.shape[0]:
        return gr[0], gr[1], np.zeros(gr.shape[1])
    else:
        return gr[0], np.zeros(gr.shape[1]), np.zeros(gr.shape[1])

def _times_to_intervals(endtimes, intervals, n):
    allpos = True
    lasttime = 0.0

    for val in range(n):
        intervals[val] = endtimes[val] - lasttime
        lasttime = endtimes[val]
        if intervals[val] <= 0:
            allpos = False 
    return allpos

def process_time_points(tp, rf_length):
    """
    THREE Cases:
		1) Single value given -> this is the interval length for all.
		2) List of intervals given.
		3) Monotonically INCREASING list of end times given.

	For all cases, the goal is for tp to have the intervals.
    """
    if type(tp) == type(0.0) or type(tp) == type(0):
        return tp * np.ones(rf_length)
    elif rf_length != tp.size:
        raise IndexError("time point length is not equal to rf length")
    else:
        ti = np.zeros(rf_length)
        if _times_to_intervals(tp, ti, rf_length):
            tp = ti
    return tp        

def process_off_resonance_arguments(df):
    if type(df) == type(0.0) or type(df) == type(0):
        return (df * np.ones(1)), 1 
    return df, df.size

def process_positions(dp):
    """
    Gets positions vectors if they exist. Zeros otherwise.
    """
    if type(dp) == type(0.0) or type(dp) == type(0):
        return dp*np.ones(1), np.zeros(1), np.zeros(1), 1
    if 3 == dp.shape[0]:
        return dp[0], dp[1], dp[2], dp[0].shape
    elif 2 == dp.shape[0]:
        return dp[0], dp[1], np.zeros(dp.shape[1]), dp[0].shape
    else:
        return dp[0], np.zeros(dp.shape[1]), np.zeros(dp.shape[1]), dp[0].size 

def process_magnetization(mx_0, my_0, mz_0, rf_length, freq_pos_count, mode):
    """
    Returns mx, my, and mz vectors allocated based on input parameters.
    """
    out_points = 1
    if 2 == mode:
        out_points = rf_length
    fn_out_points = out_points * freq_pos_count
    mx = np.zeros(fn_out_points)
    my = np.zeros(fn_out_points)
    mz = np.zeros(fn_out_points)
    if None != mx_0 and type(mx_0) != type(0.0) and type(mx_0) != type(0) and freq_pos_count == mx_0.size and freq_pos_count == my_0 and freq_pos_count == mz_0:
        for val in range(freq_pos_count):
            mx[val * out_points] = mx_0[val]
            my[val * out_points] = my_0[val]
            mz[val * out_points] = mz_0[val]
    else:
        for val in range(freq_pos_count):
            mx[val * out_points] = 0
            my[val * out_points] = 0
            mz[val * out_points] = 1
    return mx, my, mz
