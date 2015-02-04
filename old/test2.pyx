import cython
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module.
cimport numpy as np

def sum_array1(tosum):    
    return sum(tosum)

def sum_array2(tosum):
    sum = 0
    for i in xrange(tosum.shape[0]):
        sum += tosum[i]
    return sum

@cython.boundscheck(False)
def sum_array3(np.ndarray[np.float_t, ndim=1] tosum):    
    # make sure array is of correct type
    assert tosum.dtype == np.float 

    cdef int dim 
    dim = tosum.shape[0]
    cdef float sum 
    sum = 0
    for i in xrange(dim):
        sum += tosum[<unsigned int>i]
    return sum
