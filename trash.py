
import pyximport; pyximport.install()
import my_cython_functions #cython code used in this script

import rpy2 #allows me to call R functions from python
import rpy2.robjects as robj
from rpy2.robjects.numpy2ri import numpy2ri #submodules not imported automatically

import os, sys

import cPickle as pickle

import pymzml
import numpy as np
import math

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import time

def main():

    pi = robj.r['pi']
    print pi 
    print pi+2
    print pi[0]
    print pi[0]+2

    nrow = 5
    ncol = 10
    counter = 0
    binned = np.zeros((nrow, ncol), dtype="float64")
    for row in xrange(nrow):
        for col in xrange(ncol):
            binned[row, col] = counter
            counter += 1
    nele = nrow*ncol
    #rbinned = np.reshape(binned, 1)
    print binned.shape
    print binned 
    
    #vec = robj.FloatVector([1.1, 2.2, 0, 4.4, 5.5, ])
    #print numpy2ri(rbinned)
    mat = robj.r['data.frame'](numpy2ri(binned), code="ID1000")
    print mat

    

if __name__ == '__main__':
    main()
