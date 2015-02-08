import cython
import numpy as np
# "cimport" is used to import special compile-time info about the numpy module
cimport numpy as np
    
def populate_2D_binned_array(intensity_3D_full, MAX_NUM_RT_VALUES,
            rt_start, rt_end, rt_grid_size, rt_bin_width, 
            mz_start, mz_end, mz_grid_size, mz_bin_width):
    """Returns 2d array of intensity values based on specified grid

    """
    # build the 2-D array (fill with zeros to begin)
    my2Da = np.zeros((rt_grid_size, mz_grid_size), dtype=float)

    for row in intensity_3D_full:

        # for testing a smaller data set
        if MAX_NUM_RT_VALUES != None and row >= MAX_NUM_RT_VALUES: 
            break 

        # skip scans with out-of-range rt or mz values
        # include start vals but not end vals in array (avoids range prob)
        if row[0]<rt_start or row[0]>=rt_end or row[1]<mz_start or row[1]>=mz_end: 
            continue

        # rt bin num is number of steps till value is reached (round down)
        rt_bin_num = int( (row[0] - rt_start)/rt_bin_width ) 
        # if current rt is [close to] max rt then bin_num is out of range
        if rt_bin_num == rt_grid_size: 
            rt_bin_num = rt_grid_size-1 #force out of range vals into range   
        # find mz bin num in similar fashion 
        mz_bin_num = int( (row[1] - mz_start)/mz_bin_width ) 
        if mz_bin_num == mz_grid_size: 
            mz_bin_num = mz_grid_size-1 
        
        # add intensity value to bin
        my2Da[rt_bin_num, mz_bin_num] += row[2]

    return my2Da
