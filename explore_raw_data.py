
###############################################################################
## Explore mzML data (after converting .d files)
###############################################################################
import cPickle as pickle

import pymzml
import numpy as np

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import time


###############################################################################

datadir = "/media/scratch/carolyn/data_mzML/serumR1/" #on Wotan
outdir = "/home/carolyn/results_dengue/"  #on Wotan

###############################################################################
"""
class LCMSRun():
    def _find_mz_rt_bounds(self):
        pass
    
    def _build_2d_array(self):
        min_mz = self.mz_bnds[0]
        pass
    
    def __init__(self, filename):
        self._filename = filename
        self.mz_bnds, self.rt_bnds = self._find_mz_rt_bounds(self)
        self.2d_array = self._build_2d_array(self)
        self.rt_vals, self.mz_vals, self.intensity_vals = self._build_1d_arrays(self)
    
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.rt_vals, self.mz_vals, self.intensity_vals, color='r', marker='o')
        ax.set_xlabel('retention times')
        ax.set_ylabel('m/z values')
        ax.set_zlabel('Intensity (sum in bins)')
        plt.title("SerumR1: sample 1")
        plt.show() #works if you ssh -X (from the lab network, else too slow)
        plt.savefig(outdir+'scatter_mz_rt_i_'+name_suffix+'.png')
        plt.close()

    def build_spline_basis(self):
        return numpy.build_splines(self.2d_array)

def load_lcms_run(fname)
    try:
        with open(filename + ".pickled") as fp:
            return pickle.load(fp)
    except IOError:
        my_run = LCMSRun(filename)
        with open(filename + ".pickled", "w") as fp:
            pickle.dump(my_run, fp)
        return my_run

#my_run = load_lcms_run(filename)
#my_run.plot()
"""

def get_rt(ms_scan):
    for x in ms_scan.xmlTree:
        if x.get('accession') == "MS:1000016":
            return float(x.get('value'))

def plot_mz_vs_intensity(mysample):

    # obtain an iterator over scan start times
    lcms_run = pymzml.run.Reader(datadir+mysample)
    ms_scan = next(lcms_run) #data for first rt time

    mz_array = np.array(ms_scan.mz)
    i_array = np.array(ms_scan.i)
    print i_array
    plt.scatter(x=mz_array, y=i_array)
    plt.xlabel("MZ value")
    plt.ylabel("Intensity")
    plt.title("SerumR1: Sample 1")
    #plt.show() #works if you ssh -X (from the lab network, else too slow)
    #create png -- pdf is bad idea (creates large file)
    plt.savefig(outdir+'plot_mz_vs_intensity_v1.png') 
    pass

def find_mz_rt_bounds(mysample):
    # find the min and max mz values (across all retention times)
      # and the min and max retention times
    min_mz = 0.0438833333333
    max_mz = 1699.52545493
    min_rt = 0.0438833333333
    max_rt = 44.99415
    """
    # obtain an iterator over scan start times
    lcms_run = pymzml.run.Reader(datadir+mysample)
    min_mz = float(1e100)
    max_mz = 0.
    min_rt = float(1e100)
    max_rt = 0.
    for ms_scan in lcms_run:
        #Note: "None" is considered less than all numbers
        current_min_mz = min(x for x in ms_scan.mz if x!= None) 
        min_mz = min(min_mz, current_min_mz)
        max_mz = max(max_mz, max(ms_scan.mz))
        current_rt = get_rt(ms_scan)
        if current_rt != None:
            min_rt = min(min_rt, current_rt)
        max_rt = max(max_rt, current_rt)"""
    print "min mz value: " , min_mz
    print "max mz value: " , max_mz
    print "min rt value: " , min_rt #interesting -- identical to min_mz value
    print "max rt value: " , max_rt
    return min_mz, max_mz, min_rt, max_rt


def build_2D_intensity_array(
    rt_grid_size, min_rt, rt_step_size, 
    mz_grid_size, min_mz, mz_step_size, mysample):
    """Returns 2D array that contains intensity sums over grid

    """
    # build the 2-D array (fill with zeros to begin)
    mya = np.zeros((rt_grid_size, mz_grid_size), dtype=float)
    print "shape of mya" , mya.shape 

    rt_miss_count = 0
    rt_nomiss_count = 0
    mz_miss_count = 0
    i_miss_count = 0

    lcms_run = pymzml.run.Reader(datadir+mysample)
    for scan_count, ms_scan in enumerate(lcms_run):
        #if scan_count<20: #for testing  
        #rt is missing about 20% of time -- skip over these scans
        if get_rt(ms_scan) == None: 
            rt_miss_count += 1
            continue
        rt_nomiss_count += 1
        #print "rt value: " , get_rt(ms_scan)
        # bin number will be the number of steps till value is reached (round down)
        rt_bin_num = int( (get_rt(ms_scan) - min_rt)/rt_step_size ) # int = floor
        #when current rt is the max rt then rt_bin_num will be out of range
        if rt_bin_num == rt_grid_size: rt_bin_num = rt_grid_size-1 
        #ms_scan.peaks is an array of tuples
        for current_mz, current_i in ms_scan.peaks: 
            if current_mz == None:
                mz_miss_count += 1
                continue
            #print "mz value: " , current_mz
            mz_bin_num = int( (current_mz - min_mz)/mz_step_size ) #int = floor
            if mz_bin_num == mz_grid_size: mz_bin_num = mz_grid_size-1 
            #print "bin (mz, rt): " , mz_bin_num, rt_bin_num
            # increment intensity value in bin
            if current_i == None:
                i_miss_count += 1
                continue
            mya[rt_bin_num, mz_bin_num] += current_i
        #print "scan_count: " , scan_count
    print "# of non-missing rt values: " , rt_nomiss_count
    print "# of missing rt values: " , rt_miss_count
    print "# of missing mz values when rt is non-missing: " , mz_miss_count
    print "# of missing intensities when rt, mz are non-missing: " , i_miss_count
    return mya


def build_1D_arrays_from_2D_array(
            rt_grid_size, min_rt, max_rt, rt_step_size,
            mz_grid_size, min_mz, max_mz, mz_step_size, mya):
    """Returns 1D array for each of rt, mz, and intensity values

    """
    bin_count = rt_grid_size * mz_grid_size

    rt_vals = np.zeros(bin_count, dtype=float)
    bin = 0
    for row_val in np.arange(min_rt, max_rt, rt_step_size):
        for col_count in xrange(mz_grid_size):
            rt_vals[bin] = row_val
            bin += 1
    #print "rt_vals: ", rt_vals 

    mz_vals = np.array( list( np.arange(min_mz, max_mz, mz_step_size) ) * rt_grid_size )
    #print "mz_vals: ", mz_vals    

    intensity_vals = np.zeros(bin_count, dtype=float)
    bin = 0
    for row_count in xrange(rt_grid_size):
        for col_count in xrange(mz_grid_size):
            intensity_vals[bin] = mya[row_count, col_count] 
            bin += 1
    #print "intensity_vals: ", intensity_vals

    return rt_vals, mz_vals, intensity_vals
 


def build_mz_rt_intensity_arrays(mysample, mz_grid_size, rt_grid_size):
    """Returns a two-dimensional array of intensities

    """
    # find bounds
    min_mz, max_mz, min_rt, max_rt = find_mz_rt_bounds(mysample)
    mz_step_size = (max_mz - min_mz)/mz_grid_size
    print "mz_step_size: " , mz_step_size
    rt_step_size = (max_rt - min_rt)/rt_grid_size
    print "rt_step_size: " , rt_step_size

    # build 2D array with sum of intensity values within bins (reasonable??)
    intensity_2D_array = build_2D_intensity_array(
        rt_grid_size, min_rt, rt_step_size, 
        mz_grid_size, min_mz, mz_step_size, mysample)

    # obtain rt and mz values corresonding to rows and columns of 2D array
    rows_2D_array = np.arange(min_rt, max_rt, rt_step_size)
    cols_2D_array = np.arange(min_mz, max_mz, mz_step_size)

    # build 1D arrays for X, Y, Z using 2D array
    rt_vals, mz_vals, intensity_vals = build_1D_arrays_from_2D_array(
            rt_grid_size, min_rt, max_rt, rt_step_size,
            mz_grid_size, min_mz, max_mz, mz_step_size, intensity_2D_array)

    return intensity_2D_array, rows_2D_array, cols_2D_array, rt_vals, mz_vals, intensity_vals       


def plot_mz_rt_intensity(mysample, mz_grid_size, rt_grid_size, name_suffix):

    # build 2D and 1D arrays containing sum of intensities within rt/mz bins
    intensity_2D_array, rows_2D_array, cols_2D_array, rt_vals, mz_vals, intensity_vals = build_mz_rt_intensity_arrays(
        mysample, mz_grid_size, rt_grid_size)
    
    # scatter plot of the intensity array (not very clear)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(rt_vals, mz_vals, intensity_vals, color='r', marker='o')
    ax.set_xlabel('retention times')
    ax.set_ylabel('m/z values')
    ax.set_zlabel('Intensity (sum in bins)')
    plt.title("SerumR1: sample 1")
    plt.show() #works if you ssh -X (from the lab network, else too slow)
    plt.savefig(outdir+'scatter_mz_rt_i_'+name_suffix+'.png')
    plt.close()
    
    # wire frame plot (clearer than scatter)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(rt_vals, mz_vals, intensity_vals, rstride=10, cstride=10)
    ax.set_xlabel('retention times')
    ax.set_ylabel('m/z values')
    ax.set_zlabel('Intensity (sum in bins)')
    plt.title("SerumR1: sample 1")
    plt.savefig(outdir+'wireframe_mz_rt_i_'+name_suffix+'.png')
    plt.close()
    
    # heatmap
    fig, ax = plt.subplots()
    mycf = plt.contourf(rows_2D_array, cols_2D_array, intensity_2D_array)
    mycbar = plt.colorbar(mycf) #You must first define an image, eg with imshow
    mycbar.set_label("Intensities")
    ax.set_xlabel('Retention times')
    ax.set_ylabel('m/z values')
    plt.title("SerumR1: sample 1")
    plt.savefig(outdir+'contour_mz_rt_i_'+name_suffix+'.png')    
    plt.close()

    return


###############################################################################

def main():

    start_time = time.time()
    
    mysample = "Nicaserhilic1000DF5d.mzML" 
   
    # plot mz vs. intensity for just 1 scan start time
    #plot_mz_vs_intensity(mysample) 

    # plot mz vs. rt vs. intensity for various grid resolutions 
    plot_mz_rt_intensity(
        mysample, mz_grid_size=20, rt_grid_size=20, name_suffix="20by20")
    plot_mz_rt_intensity(
        mysample, mz_grid_size=100, rt_grid_size=100, name_suffix="100by100")  
    
    print "Execution time: " , time.time() - start_time, "seconds"

if __name__ == '__main__':
    main()


###############################################################################
################################## Scratch code ###############################
###############################################################################

def scratch():
    for x in t.xmlTree: print x.get('accession'), x.get('name'), x.get('value')

    for ms_scan in lcms_run:
        illution_time = get_rt(ms_scan) #gives just one rt value
        mz_array = ms_scan.mz #gives an entire array of mz values
        intensity_array = ms_scan.i #gives an entire array of intensities
        print illution_time
        print mz_array
        print intensity_array
    pass


