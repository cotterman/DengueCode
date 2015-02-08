
import pyximport; pyximport.install()
import my_cython_functions #cython code used in this script

import rpy2 #allows me to call R functions from python
import rpy2.robjects as robj
from rpy2.robjects.numpy2ri import numpy2ri #submodules not imported automatically
#the following import will allow me to import arbitrary R code as a package
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

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

def screw_around():

    pi = robj.r['pi']
    print pi 
    print pi+2
    print pi[0]
    print pi[0]+2

    #create fake binned array
    nrow = 5
    ncol = 10
    counter = 0
    binned = np.zeros((nrow, ncol), dtype="float64")
    for row in xrange(nrow):
        for col in xrange(ncol):
            binned[row, col] = counter
            counter += 1
    #print binned
    
    #get binned array into R data.frame
    #vec = robj.FloatVector([1.1, 2.2, 0, 4.4, 5.5, ])
    #print binned.shape
    print numpy2ri(binned)
    rdf = robj.r['data.frame'](numpy2ri(binned), code="ID1000")
    #print rdf

    # now see if we can get R to use this dataframe 
    myRcode = """
    square <- function(rdf) {
        myv = rdf$X2 + rdf$X3
        return(myv)
    }
    doit <- function() {
        source("/srv/scratch/carolyn/Dengue_code/Rtest_rpy.R") 
        run_test_wrap(3)
    }
    """
    print "wwwwah"
    powerpack = SignatureTranslatedAnonymousPackage(myRcode, "powerpack")
    print powerpack._rpy2r.keys() #to reveal the functions within powerpack
    print powerpack.square(rdf) #to run the function "square" found in powerpack
    print powerpack.doit()


def build_row_of_lcms_matrix(binned, respD, nrow, ncol, filecount, filename):

    print os.path.abspath(filename)
    # First column should contain "code", which is IDXXXX
    start_pos = filename.find("Nicaserhilic")
    IDcode = "ID" + filename[start_pos+12:start_pos+16]
    print "ID: " , IDcode
    respD[filecount,0] = IDcode
    cell = 1
    for row in xrange(nrow):
        for col in xrange(ncol):
            respD[filecount,cell] = binned[row,col]
            cell += 1
    return respD
    
    
###############################################################################

def parse_arguments():
    """ When running from command prompt, expect filename and output directory

    Ex: python /srv/scratch/carolyn/Dengue_code/process_raw_data_and_do_prediction.py
               /srv/scratch/carolyn/mzml_serumR1/Ni*.mzML
               /srv/scratch/carolyn/Results/
    """ 
    return sys.argv[1:-1], os.path.abspath(sys.argv[-1])

def main():
    filenames, outdir = parse_arguments() #filenames will be a list
    os.chdir(outdir) #change pwd to output directory
    #print "filenames", filenames

    # will need to get intensity_2D_binned into R data.frame
        # intensity_2D_binned was created in my_cython_functions.pyx as follows:
            # cdef np.ndarray[np.float_t, ndim=2] my2Da
            # my2Da = np.zeros((rt_grid_size, mz_grid_size)) """
    #create fake binned array
    nrow = 5
    ncol = 10
    counter = 0
    binned = np.zeros((nrow, ncol), dtype="float64")
    for row in xrange(nrow):
        for col in xrange(ncol):
            binned[row, col] = counter
            counter += 1

    ### Option 1 ###
    #turn each binned array into one row of what will be an R data.frame
        # then add ID to each row and combine these binned arrays
    nele = nrow*ncol
    rbinned1 = np.reshape(binned, nele)
    rbinned1c = np.hstack(( np.array(["ID1001"]), rbinned1 )) #concatenate also works
    rbinned2 = np.reshape(binned, nele)
    rbinned2c = np.hstack(( np.array(["ID1002"]), rbinned2 ))
    stacked = np.vstack((rbinned1c, rbinned2c))
    #print stacked
    #print stacked.shape #2 by 51
    rrdf1 = robj.r['data.frame'](numpy2ri(stacked))
    #print rrdf1

    ### Option 2 ###
    #build empty array that is 5 (# mzml files) by 51 (# rt/mz bins + 1 for patient ID)
        # fill each row with the binned data
    floatD = np.zeros((5,nrow*ncol), dtype="float64")
    strD = np.zeros((5,1), dtype='a6') #a6 is the dtype for a 6 character string
    respD = np.hstack((strD, floatD))
    print respD.shape

    for filecount, filename in enumerate(filenames):
        if filecount<2:
            respD = build_row_of_lcms_matrix(
                binned, respD, nrow, ncol, filecount, filename)
    #print respD
    df2 = robj.r['data.frame'](numpy2ri(respD))
    print df2

    # now see if we can get R to use this dataframe 
    myRcode = """
    doR <- function(python_respD, lcms_run) {
        source("/srv/scratch/carolyn/Dengue_code/prediction_with_LCMS_from_python.R") 
        run_predictions_wrap(python_respD, lcms_run)
    }
    """
    #Rpack = SignatureTranslatedAnonymousPackage(myRcode, "Rpack")
    #print Rpack._rpy2r.keys() #to reveal the functions within powerpack
    #3print Rpack.doR(df2, 1) #to run the function found in powerpack

    # now see if we can get R to use this dataframe 
    myRcode = """
    square <- function(rdf) {
        myv = rdf$X2 + rdf$X3
        return(myv)
    }
    doit <- function(input) {
        source("/srv/scratch/carolyn/Dengue_code/Rtest_rpy.R") 
        run_test_wrap(input)
    }
    """
    print "wwwwah"
    powerpack = SignatureTranslatedAnonymousPackage(myRcode, "powerpack")
    #print powerpack._rpy2r.keys() #to reveal the functions within powerpack
    #print powerpack.square(df2) #to run the function "square" found in powerpack
    print powerpack.doit(df2)

   

if __name__ == '__main__':
    main()

