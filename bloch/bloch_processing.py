import numpy as np
import scipy as sp

NUMBER = (int, float, complex)

#Functions to handle preprocessing for bloch simulator arguments.

def process_gradient_argument(gr, points):
    """
    Takes in a gradient argument and returns directional gradients.
    If gradients don't exist, returns array of zeros.
    If only one number is passed, it's assigned to the entire xgrad and the others are zeros.
    """
    if isinstance(gr, NUMBER):
        return gr * np.ones(points), np.zeros(points), np.zeros(points)
    elif 1 == len(gr.shape):
        return gr, np.zeros(points), np.zeros(points)

    gradient_dimensions = gr.shape[0]

    if 3 == gradient_dimensions:
        return gr[0,:], gr[1,:], gr[2,:]
    elif 2 == gradient_dimensions:
        return gr[0,:], gr[1,:], np.zeros(points)
    else:
        return gr[0,:], np.zeros(points), np.zeros(points)

def process_time_points(tp, points):
    """
    THREE Cases:
		1) Single value given -> this is the interval length for all.
		2) List of intervals given.
		3) Monotonically INCREASING list of end times given.

	For all cases, the goal is for tp to have the intervals.
    """
    if isinstance(tp, NUMBER):
        return tp * np.ones(points)
    elif points != tp.size:
        raise IndexError("time point length is not equal to rf length")
    else:
        ti = np.zeros(points)
        if _times_to_intervals(tp, ti, points):
            tp = ti
    return tp        

def process_off_resonance_arguments(df):
    """
    Processes off resonance arguments.
    Returns df and size. If only one numer is passed, returns number as single array.
    """
    if isinstance(df, NUMBER):
        return (df * np.ones(1)), 1 
    return df, df.size

def process_positions(dp):
    """
    Gets positions vectors if they exist. Zeros otherwise.
    If only one number is passed, is set as xgrad and other directions are 0s.
    """
    if isinstance(dp, NUMBER):
        return dp*np.ones(1), np.zeros(1), np.zeros(1), 1
    elif 1 == len(dp.shape):
        return dp, np.zeros(dp.size), np.zeros(dp.size), dp.size

    position_dimensions = dp.shape[0]
    number_of_positions = dp.shape[1]

    if 3 == position_dimensions:
        return dp[0,:], dp[1,:], dp[2,:], number_of_positions
    elif 2 == position_dimensions:
        return dp[0,:], dp[1,:], np.zeros(number_of_positions), number_of_positions
    else:
        return dp[0,:], np.zeros(number_of_positions), np.zeros(number_of_positions), number_of_positions 

def process_magnetization(mx_0, my_0, mz_0, rf_length, freq_pos_count, mode):
    """
    Returns mx, my, and mz vectors allocated based on input parameters.
    """
    if isinstance(mx_0, np.ndarray) and isinstance(my_0, np.ndarray) and isinstance(mz_0, np.ndarray):
        mx_0 = mx_0.ravel()
        my_0 = my_0.ravel()
        mz_0 = mz_0.ravel()
    out_points = 1
    if (2 & mode):
        out_points = rf_length
    fn_out_points = out_points * freq_pos_count
    mx = np.zeros(fn_out_points)
    my = np.zeros(fn_out_points)
    mz = np.zeros(fn_out_points)
    if None != mx_0 and type(mx_0) != type(0.0) and type(mx_0) != type(0) and freq_pos_count == mx_0.size and freq_pos_count == my_0.size and freq_pos_count == mz_0.size:
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

def reshape_matrices(mx, my, mz, ntime, n_pos, nf):
    """
    Reshapes output matrices.
    """
    if ntime > 1 and nf > 1 and n_pos > 1:
        shape = (nf, n_pos, ntime)
        mx.shape = shape
        my.shape = shape
        mz.shape = shape
        return
    else:
        if ntime > 1:
            shape = ((n_pos * nf), ntime)
            if  1 == (n_pos * nf):
                shape = (ntime, )
        else:
            shape = (nf, n_pos)
            if 1 == nf:
                shape = (n_pos,)
        mx.shape = shape
        my.shape = shape
        mz.shape = shape

def _times_to_intervals(endtimes, intervals, n):
    """
    Helper function for processing time points.
    """
    allpos = True
    lasttime = 0.0

    for val in range(n):
        intervals[val] = endtimes[val] - lasttime
        lasttime = endtimes[val]
        if intervals[val] <= 0:
            allpos = False 
    return allpos


