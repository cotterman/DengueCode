import cython
import numpy as np
# "cimport" is used to import special compile-time info about the numpy module
cimport numpy as np

@cython.boundscheck(False)  
@cython.cdivision(True)  
def populate_2D_binned_array(
            np.ndarray[np.float_t, ndim=2] intensity_3D_full, 
            MAX_NUM_RT_VALUES,
            np.float_t rt_start, np.float_t rt_end, np.int_t rt_grid_size, 
            np.float_t mz_start, np.float_t mz_end, np.int_t mz_grid_size):
    """Returns 2d array of intensity values based on specified grid

    """
    assert rt_grid_size > 0 and mz_grid_size > 0
    assert rt_end > rt_start and mz_end > mz_start

    cdef unsigned int rt_bin_num
    cdef unsigned int mz_bin_num
    cdef unsigned int c_MAX_NUM_RT_VALUES

    # calculate bin widths
    cdef np.float_t rt_bin_width
    cdef np.float_t mz_bin_width 
    rt_bin_width = (rt_end - rt_start)/rt_grid_size 
    mz_bin_width = (mz_end - mz_start)/mz_grid_size

    # build the 2-D array (fill with zeros to begin)
    cdef np.ndarray[np.float_t, ndim=2] my2Da
    my2Da = np.zeros((rt_grid_size, mz_grid_size))

    # in C, it doesn't make sense for an integer to equal None
    if MAX_NUM_RT_VALUES == None:
        c_MAX_NUM_RT_VALUES = intensity_3D_full.shape[0] + 1
    else:
        c_MAX_NUM_RT_VALUES = MAX_NUM_RT_VALUES

    cdef unsigned int row_num
    for row_num in xrange(intensity_3D_full.shape[0]):

        # for testing a smaller data set
        if row_num >= c_MAX_NUM_RT_VALUES: 
            break 

        # skip scans with out-of-range rt or mz values
        # include start vals but not end vals in array (avoids range prob)
        if (intensity_3D_full[row_num,0] < rt_start 
             or intensity_3D_full[row_num,0] >= rt_end
             or intensity_3D_full[row_num,1] < mz_start 
             or intensity_3D_full[row_num,1] >= mz_end): 
            continue

        # rt bin num is number of steps till value is reached (round down)
        rt_bin_num =  <unsigned int>(
            (intensity_3D_full[row_num,0] - rt_start)/rt_bin_width)
        # if current rt is [close to] max rt then bin_num is out of range
        if rt_bin_num == rt_grid_size: 
            rt_bin_num = rt_grid_size-1 #force out of range vals into range   

        # find mz bin num in similar fashion 
        mz_bin_num = <unsigned int>(
            (intensity_3D_full[row_num,1] - mz_start)/mz_bin_width) 
        if mz_bin_num == mz_grid_size: 
            mz_bin_num = mz_grid_size-1 
        
        # add intensity value to bin
        my2Da[rt_bin_num, mz_bin_num] += intensity_3D_full[row_num,2]

    return my2Da
