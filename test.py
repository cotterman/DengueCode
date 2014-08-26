
# This line tells python to compile the .pyx code before running
    # to compile the .pyx code and not run the whole python script:
    # $ cython -a my_cython_code.pyx
import pyximport; pyximport.install()
import test2
import numpy as np


tosum = np.random.randn(10000000)
print type(tosum)
test2.sum_array3(tosum)
#tosum.sum()
