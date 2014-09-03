
###############################################################################
## Process mzML data -- build arrays and graphs (after converting .d files)  ##
###############################################################################

# This line tells python to compile the .pyx code before running
    # to compile the .pyx code and not run the whole python script:
    # $ cython -a my_cython_code.pyx
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

# limit the number of data points that we read, mostly for debugging
MAX_NUM_RT_VALUES = 500 #set to None to run on all data

VERBOSE = True
def log_statement(statement):
    if VERBOSE: print statement

###############################################################################

def load_saved_LCMSRun(fname):
    with open(fname) as fp:
        lcms_run = pickle.load(fp)
    with open(fname + ".mat") as fp:
        mat = np.load(fp)
    lcms_run.intensity_3D_full = mat
    return lcms_run
    
class LCMSRun():
    'Class for mzML files (raw LCMS data)'

    def save(self, ofname):
        #save the full 3D intensity array into a .mat file
            # we're doing this b/c the pickled 3D array was huge and slow 
        with open(ofname + ".mat", "w") as ofp:
            np.save(ofp, self.intensity_3D_full)
        intensity_3D_full = self.intensity_3D_full
        self.intensity_3D_full = None
        #dump the rest of the object into the pickle file
        with open(ofname, "w") as ofp:
            pickle.dump(self, ofp)
        self.intensity_3D_full = intensity_3D_full

    @staticmethod
    def get_rt(ms_scan):
        """Takes one item in an LCMSRun iterator and returns retention time

        """
        for x in ms_scan.xmlTree:
            if x.get('accession') == "MS:1000016":
                return float(x.get('value'))
    
    @staticmethod
    def get_min_mz(ms_scan):
        """Takes one item in an LCMSRun iterator and returns min mz val

        Not using this function b/c it sometimes returns None
        """
        for x in ms_scan.xmlTree:
            if x.get('accession') == "MS:1000528":
                return float(x.get('value'))

    @staticmethod    
    def get_max_mz(ms_scan):
        """Takes one item in an LCMSRun iterator and returns max mz val

        """
        for x in ms_scan.xmlTree:
            if x.get('accession') == "MS:1000527":
                return float(x.get('value'))

    @staticmethod
    def plot_mz_vs_intensity(ms_scan):
        """Takes one item in an LCMSRun iterator and plots mz v intensity

        """
        mz_array = np.array(ms_scan.mz)
        i_array = np.array(ms_scan.i)
        plt.scatter(x=mz_array, y=i_array)
        plt.xlabel("MZ value")
        plt.ylabel("Intensity")
        plt.title("MZ vs. Intensity for one retention time value")
        plt.savefig(outdir+'plot_mz_vs_intensity.png') 
    
    def _find_mz_rt_bounds(self):
        # find the min and max mz values (across all retention times)
          # and the min and max retention times
        """
        min_mz = 0.0438833333333
        max_mz = 1699.52545493
        min_rt = 0.0438833333333
        max_rt = 44.99415
        """
        # obtain an iterator over scan start times
        lcms_run = pymzml.run.Reader(self.filepath)
        min_rt = float(1e100)
        max_rt = 0.
        min_mz = float(1e100)
        max_mz = 0.
        for ms_scan in lcms_run:
            #Note: "None" is considered less than all numbers
            current_min_mz = min(x for x in ms_scan.mz if x!= None) 
            min_mz = min(min_mz, current_min_mz)
            max_mz = max(max_mz, self.get_max_mz(ms_scan) )
            current_rt = self.get_rt(ms_scan)
            if current_rt != None:
                min_rt = min(min_rt, current_rt)
            max_rt = max(max_rt, current_rt)
        print "min mz value: " , min_mz
        print "max mz value: " , max_mz
        print "min rt value: " , min_rt # identical to min_mz value
        print "max rt value: " , max_rt
        #return 2 tuples (note: tuples are immutable)
        
        return (min_rt, max_rt), (min_mz, max_mz)

    def _find_num_data_points(self):
        lcms_run = pymzml.run.Reader(self.filepath)
        num = 0
        for scan_num, ms_scan in enumerate(lcms_run):
            if MAX_NUM_RT_VALUES != None and scan_num >= MAX_NUM_RT_VALUES: 
                break
            num += ms_scan.get('defaultArrayLength')
        
        return num
    
    def _build_3D_full_array(self):
        """Returns 3D array of intensity values without summarization

            Also returns array of rt,mz,i values when at least 1 is missing
        """
        # build the 3-D array (fill with zeros to begin)
        # the size of the array will need to equal the number of rt/mz combos
        num_rt_mz_combos = self._find_num_data_points()
        log_statement("Number of rt/mz combos: {}".format(num_rt_mz_combos))
        
        # allocate space to hold all mz/rt/intensity values
        data = np.zeros((num_rt_mz_combos, 3), dtype=float)
        missing_values = []
        
        # this array will have missing values where LCMS data has missings
        log_statement("Populating full array")
        lcms_run = pymzml.run.Reader(self.filepath)
        row = 0
        for ms_count, ms_scan in enumerate(lcms_run):
            # for testing a smaller data set
            if MAX_NUM_RT_VALUES != None and ms_count >= MAX_NUM_RT_VALUES: 
                break 
            
            rt = self.get_rt(ms_scan)
            for mz, intensity in ms_scan.peaks:
                try:
                    rt, mz, intensity = float(rt), float(mz), float(intensity)
                except TypeError:
                    missing_values.append((rt, mz, intensity))
                    #print ms_count, row, missing_values
                    continue
                data[row, 0] = rt
                data[row, 1] = mz
                data[row, 2] = math.log(intensity+1)
                row += 1

        # resize the data array because some observations may have been skipped 
        # due to missing values
        data.resize((row, 3))
        return data, missing_values

    # Create a class constructor/initiation method 
    # (This is called when you create a new instance of this class)
    def __init__(self, filepath):
        self.filepath = filepath
        log_statement("Finding mz and rt bounds")
        # can we speed this up by looking through the XML file?
        self.rt_bounds, self.mz_bounds = self._find_mz_rt_bounds()
        log_statement("Building 3D array")
        ( self.intensity_3D_full, self.missing_values 
          ) = self._build_3D_full_array()
        print "Missing values: " , self.missing_values

    def build_2D_binned_array(self, MAX_NUM_RT_VALUES, 
            rt_start, rt_end, rt_grid_size, 
            mz_start, mz_end, mz_grid_size):
        """Returns 2d array of intensity values based on specified grid

        """

        # populate 2D binned array using full intensity array
        return my_cython_functions.populate_2D_binned_array(self.intensity_3D_full, 
            MAX_NUM_RT_VALUES,
            rt_start, rt_end, rt_grid_size,  
            mz_start, mz_end, mz_grid_size)

    def make_marginal_histograms(self):
        """ Create histograms of rt, mz, and intensity values in data

        """
        """
        fig, ax = plt.subplots()
        myhist_rt = plt.hist(self.intensity_3D_full[:,0])
        ofname = ('hist_rt.png')
        plt.savefig(ofname)
        
        fig, ax = plt.subplots()
        myhist_mz = plt.hist(self.intensity_3D_full[:,1])
        ofname = ('hist_mz.png')
        plt.savefig(ofname)
        """
        fig, ax = plt.subplots()
        """
        #only plot for intensities greater than 200,000 (takes 7 min.)
        high_intensities = np.zeros(self.intensity_3D_full.shape[0])
        counter = 0
        for intensity in self.intensity_3D_full[:,2]:
            if intensity > 200000:
                high_intensities[counter] = intensity
                counter +=1
        high_intensities.resize(counter)
        myhist_i = plt.hist(high_intensities)
        """
        myhist_i = plt.hist(self.intensity_3D_full[:,2])
        ofname = ('hist_log_i.png')
        plt.savefig(ofname)
        
        return

def build_1D_binned_arrays(intensity_2D_binned, 
        rt_start, rt_end, rt_grid_size, mz_start, mz_end, mz_grid_size):
    """Returns 1D array for each of rt, mz, and i vals using 2D array

    """

    bin_count = rt_grid_size * mz_grid_size

    #the bin value will correspond to the bin's start value
    rt_bin_val = np.linspace(rt_start, rt_end, rt_grid_size, endpoint=False)
    rt_1D_binned = np.zeros(bin_count, dtype=float)
    bin = 0
    for row_val in rt_bin_val:
        for col_count in xrange(mz_grid_size):
            rt_1D_binned[bin] = row_val
            bin += 1
    #print "rt_1D_binned: ", rt_1D_binned 
    print "# of rt_1D_binned: " , rt_1D_binned.shape

    #the mz bin value will correspond to the bin's start value
    mz_bin_val = np.linspace(mz_start, mz_end, mz_grid_size, endpoint=False)
    mz_1D_binned = np.array( list( mz_bin_val ) * rt_grid_size )
    #print "mz_1D_binned: ", mz_1D_binned    
    print "# of mz_1D_binned: " , mz_1D_binned.shape

    intensity_1D_binned = np.zeros(bin_count, dtype=float)
    bin = 0
    for row_count in xrange(rt_grid_size):
        for col_count in xrange(mz_grid_size):
            intensity_1D_binned[bin] = (
                intensity_2D_binned[row_count, col_count]) 
            bin += 1
    #print "intensity_1D_binned: ", intensity_1D_binned
    print "# of intensity_1D_binned: " , intensity_1D_binned.shape

    return rt_1D_binned, mz_1D_binned, intensity_1D_binned
    

def scatter3D(ofname, rt_1D_binned, mz_1D_binned, intensity_1D_binned):
    """Generate scatter plot using the 1D intensity arrays 

    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(rt_1D_binned, mz_1D_binned, intensity_1D_binned, 
        color='r', marker='o')
    ax.set_xlabel('retention times')
    ax.set_ylabel('m/z values')
    ax.set_zlabel('Intensity (sum in bins)')
    plt.title("SerumR1: Sample 2")
    plt.savefig(ofname)
    plt.close()
    return

def wireframe3D(ofname, rt_1D_binned, mz_1D_binned, intensity_1D_binned):
    """Generate wireframe plot using the 1D intensity arrays 

    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(rt_1D_binned, mz_1D_binned, intensity_1D_binned, 
        rstride=10, cstride=10)
    ax.set_xlabel('retention times')
    ax.set_ylabel('m/z values')
    ax.set_zlabel('Intensity (sum in bins)')
    plt.title("SerumR1: Sample 2")
    plt.savefig(ofname)
    plt.close()
    return

def heatmap(ofname, rt_bin_val, mz_bin_val, intensity_2D_binned):
    """Generate heatmap using the 2D intensity array 

    """
    fig, ax = plt.subplots()
    mycf = plt.contourf(rt_bin_val, mz_bin_val, intensity_2D_binned)
    mycbar = plt.colorbar(mycf) 
    mycbar.set_label("Intensities")
    ax.set_xlabel('Retention times')
    ax.set_ylabel('m/z values')
    plt.title("SerumR1: Sample 2")
    plt.savefig(ofname)    
    plt.close()
    return


def plot_binned_data(intensity_2D_binned, 
        rt_start, rt_end, rt_grid_size, mz_start, mz_end, mz_grid_size):
    """Get 1D arrays of binned rt and mz values and then plot binned data

    """
    # the scatter3D and wireframe3D functions require these 1D arrays as input
    log_statement("creating 1D binned arrays") 
    checkpt_start1Dbin = time.time()
    rt_1D_binned, mz_1D_binned, intensity_1D_binned = (
        build_1D_binned_arrays(intensity_2D_binned, 
            rt_start, rt_end, rt_grid_size, mz_start, mz_end, mz_grid_size) )
    rt_bin_val = np.linspace(rt_start, rt_end, rt_grid_size, endpoint=False)
    mz_bin_val = np.linspace(mz_start, mz_end, mz_grid_size, endpoint=False)
    checkpt_1Dbin = time.time()
    log_statement( "Time to create 1D binned arrays: {} minutes".format(
        (checkpt_1Dbin - checkpt_start1Dbin)/60. ) ) #negligible time
    
    log_statement("plotting the 3D data")
    scatter3D("scatter.png", rt_1D_binned, mz_1D_binned, intensity_1D_binned)
    wireframe3D("wireframe.png", rt_1D_binned, mz_1D_binned, intensity_1D_binned)
    heatmap("heatmap.png", rt_bin_val, mz_bin_val, intensity_2D_binned)
    checkpt_3Dplots = time.time()
    log_statement("Time to create 3D graphs: {} minutes".format(
        (checkpt_3Dplots - checkpt_1Dbin)/60. ) ) #takes only a second

    return

def load_lcms_run(filepath, replace=0):
    """Create or load class lcms_run instance

    """
    pickled_fname = filepath + ".pickled"
    # If pickled class instance already exists, then use it
    if replace == False:
        try:
            return load_saved_LCMSRun(pickled_fname)
        except IOError:
            pass
    
    my_run = LCMSRun(filepath) 
    my_run.save(pickled_fname)
    return my_run


def fill_row_of_lcms_matrix(respD, rt_grid_size, mz_grid_size, filecount, filename):

    filepath = os.path.abspath(filename)
    start_time = time.time()

    # Create full 3D data array and get basic stats as part of class attributes
    my_run = load_lcms_run(filepath, replace=False)
    checkpt_load = time.time()
    log_statement( "Time to build/load full 3D data array: {} minutes".format(
            (checkpt_load - start_time)/60. ) )  
    # takes ~14 min to build full array; ~20 sec to load
    log_statement("rt bounds: {}".format(my_run.rt_bounds)) # no () for attributes
    log_statement("mz bounds: {}".format(my_run.mz_bounds))
    
    #Get histograms of all the intensity, mz, and rt values in data
    #my_run.make_marginal_histograms()
    checkpt_hist = time.time()
    #log_statement("Time to create histograms: {} minutes".format(
    #    (checkpt_hist - checkpt_load)/60. ) ) #~30 seconds
   
    # Create binned data arrays with specified range and grid sizes
        # This array will be of dimension rt_grid_size by mz_grid_size
    log_statement("creating 2D binned array") 
    rt_start = 0    #rt_start = my_run.rt_bounds[0]
    rt_end = 45     #rt_end = my_run.rt_bounds[1]
    mz_start = 100  #mz_start = my_run.mz_bounds[0]
    mz_end = 1700   #mz_end = my_run.mz_bounds[1]
    intensity_2D_binned = my_run.build_2D_binned_array(MAX_NUM_RT_VALUES,
        rt_start, rt_end, rt_grid_size, mz_start, mz_end, mz_grid_size)
    checkpt_2Dbin = time.time()
    # used to take ~15 min. with pure python; now with cypython takes ~20 sec 
    log_statement( "Time to create 2D binned array: {} minutes".format(
        (checkpt_2Dbin - checkpt_hist)/60. ) ) 

    # Plots of  binned data (mz vs. rt vs. intensity)
    #plot_binned_data(intensity_2D_binned, 
    #    rt_start, rt_end, rt_grid_size, mz_start, mz_end, mz_grid_size)
    
    # First column should contain "code", which is IDXXXX
    start_pos = filename.find("Nicaserhilic")
    IDcode = "ID" + filename[start_pos+12:start_pos+16]
    print "ID: " , IDcode
    respD[filecount,0] = IDcode
    cell = 1
    for row in xrange(rt_grid_size):
        for col in xrange(mz_grid_size):
            respD[filecount,cell] = intensity_2D_binned[row, col]
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
    start_time_overall = time.time()
  
    rt_grid_size = 50
    mz_grid_size = 50
    log_statement( "Number of mzml patient files: {}".format(len(filenames)) )

    #build empty np array to be filled with LCMS values for R prediction
        #dimension will be (# mzml files) by (# rt/mz bins + 1 for patient ID)
    floatD = np.zeros((len(filenames),rt_grid_size*mz_grid_size), dtype=float) #intensity vals
    strD = np.zeros((len(filenames),1), dtype='a6') #a6 is dtype for 6 char str
    respD = np.hstack((strD, floatD))

    #fill the array
    for filecount, filename in enumerate(filenames):
        if filecount<10:
            respD = fill_row_of_lcms_matrix(respD, rt_grid_size, mz_grid_size, filecount, filename)
    print "Data for use in R: " , respD[:,0:5]

    #convert numpy array into data.frame recognized by R
    Rdf = robj.r['data.frame'](numpy2ri(respD))
    #use this dataframe in R prediction code 
    myRcode = """
    doR <- function(python_respD, lcms_run) {
        source("/srv/scratch/carolyn/Dengue_code/prediction_with_LCMS_from_python.R") 
        run_predictions_wrap(python_respD, lcms_run)
    }
    """
    Rpack = SignatureTranslatedAnonymousPackage(myRcode, "Rpack")
    print Rpack.doR(Rdf, 1) #to run the function doR, found in Rpack
      
    log_statement("Total execution time: {} minutes".format(
        (time.time() - start_time_overall)/60. ) )

if __name__ == '__main__':
    main()

# Now you can run from python:
    # import pickle
    # from process_raw_data import LCMSRun
    # datadir = "/media/scratch/carolyn/data_mzML/serumR1/"
    # my_pickle = open(datadir+'Nicaserhilic1000DF5d.mzML.pickled', 'r') #r is for read
    # my_run = pickle.load(my_pickle)




