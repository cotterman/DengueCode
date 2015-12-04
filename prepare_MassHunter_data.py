###############################################################################
####  Prepare LCMS mass hunter data for use in prediction_in_python.py ########
###############################################################################

import os, sys
import re
import cPickle as pickle
import numpy as np
import scipy as sp
import pandas as pd


def parse_arguments():
    """ When running from command prompt, expect filename and output directories

    Ex: python /users/ccotter/git_dengue/prepare_MassHunter_data.py
               /srv/scratch/ccotter/intermediate_data/
               /srv/scratch/ccotter/py_out/
    """ 
    return sys.argv[1], os.path.abspath(sys.argv[-1])


def main():
    
    inDir, outDir = parse_arguments() #filenames will be a list
    inFname = "comboD1_filter50n_wImpRF1"
    outFname = "NP_MassHunt" #name of numpy array to create
        
    # prepare data so it will properly merge with clinical data on code, Study, Cod_Nin
    NP_MH = pd.read_csv(inDir + inFname + '.txt', sep='\t')
    #print "Column names: " , list(NP_MH)
    print "Number of obs and of vars: " , NP_MH.shape
    # keep only LCMS variables and vars for merging
    keepervars = ["code","Study"]
    for var in NP_MH.columns.values:
        y = re.findall("MZ_.*", var)
        assert len(y) in (0,1)
        if len(y)==1: keepervars.append(y[0])
    print "Number of MZ values: " , len(keepervars) - 2
    NP_MH = NP_MH[keepervars]
    NP_MH.insert(2, 'Cod_Nin', np.nan) #need Cod_Nin since it will be expected later on   

    # save for later use (employ pd.read_pickle() to load)
    print "Creating pickled data frame"
    NP_MH.to_pickle(os.path.join(outDir, outFname)) 


if __name__ == '__main__':
    main()
