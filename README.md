
###############################################################################
##################### Denge Dx Project Documentation ##########################
###############################################################################

Locations
=========

## Project code: Bitbucket repository
* https://bitbucket.org/cotterman/dengue-dx 
    * this one is to be retired, as bitbucket is improper location for non-code items
* https://bitbucket.org/cotterman/dengueDx-code
    * this one is active. As name implies, contains only code-related files
    * synced to the following folders:
        * on my PC: /home/carolyn/Dengue_code (this is the E drive from Window's perspective)
        * on Amold: /srv/scratch/carolyn/dengue_mount/Dengue_code

Current organization (in which I have git code on my PC and also on Amold)
may cause confusion b/c I need to be sure to have only one working version at a time -- 
must access Wotan's mounted drive and do a git commit/push/pull to get it onto my PC afterwards.


## Location of raw data (and other large dengue data files) on server: 
* wotan.lbl.gov: 
    * /srv/scratch/carolyn/dengue_project/raw_data (contains all raw .d files)
* amold.lbl.gov:
    * /srv/scratch/carolyn/converted_serumR1, ../converted_serumR2 (contains serum mzML files)

* Note: The datasets I used for the R33 analysis are what Natalia put in folders labeled:
    * incoming/First batch 2012 (I relabled "/serumR1")
    * incoming/Second batch 2013 (I relabeled "/serumBadR2")
        * this data is bad and should probably not be used
    * incoming/Third batch 2014 (I relabeled "/serumR2)

## Location of other project materials:
* on my PC: /home/carolyn/Dengue_dx (E drive from Window's perspective)
* on Amold: /srv/scratch/carolyn


Step 1: convert .d files to mzML format 
========================================

Note: msconvert works for many conversions, but for .d, I needed the Agilent DLLS, which only works on Windows, so these 2 lines did not work:
  $ cd /media/scratch/carolyn/incoming/serumR1
  $ ~/bin/msconvert Nicaserhilic1000DF5d.d -o /media/scratch/carolyn/data_mzML/serumR1

## Transfer to PC with windows since Agilent DLLs only work on Windows:
$ cd ~/dengue_data_transfers/agilent_raw_format
$ scp -r carolyn@wotan.lbl.gov:/media/scratch/carolyn/incoming/serumR1/Nicaserhilic1000DF5d.d ./

## From windows use MSConvert (which is a ProteoWizard tool):
type CMD from the search bar (from start menu).
$ cd C:\Program Files\ProteoWizard\ProteoWizard 3.0.6337
$ msconvert --zlib --gzip --mzML C:\Users\carolyn\Documents\dengue_temp\Nicaserhilic1000DF5d.d -o C:\Users\carolyn\Documents\dengue_temp
[Conversion takes about 30 minutes per .d file]

## Transfer from PC to Server
$ cd ~/dengue_dx/Dengue_data/Raw/mzML_raw
$ scp Nicaserhilic1000DF5d.mzML.gz carolyn@wotan.lbl.gov:/media/scratch/carolyn/data_mzML/serumR1

## Unzip mzML file on server
gunzip Nicaserhilic1000DF5d.mzML.gz

## Delete files on my PC



Step 2: Create features and run prediction analysis
===================================================

## Option 1: use MZmine2 software

### Can use MZmine2's GUI or batch mode:

* Using MZmine2 GUI
    ```
    cd /dengue_dx/MZmine-2.10
    sh startMZmine_Linux.sh
    ```
    * Import raw data (takes about 10 minutes to load one mzML file)
    * Visualization: view TIC/XIC visualizer (gives rt vs. intensity) 
       * RT range: 0 - 45
       * mz range: 100 - 1800

* Using MZmine2's batch mode
    ssh -X amold

    * Create batch queue file
       * From GUI, select "batch mode" under the Project menu
       * Add the steps you with to implement
       * Click on "save" to create a .xml file with commands
    * Run this saved .xml from command line as follows
    ```
    cd ~/dengue_dx/MZmine-2.10
    sh startMZmine_Linux.sh ~/batch_queue_file.xml
    ```
    Export XML file containing resulting peak list:
    * This file will contain a <row> element for each peak
    * Each of these rows contains information on its individual peaks

### Details on MZmine2 data processing options

* scan-by-scan filtering (not necessary for high resolution data)
   * Choose mean filter with m/z bin size no less than .01
* **Baseline correction (not necessary for high resolution data)
   * chromatogram type: base peak intensity
   * smoothing factor: 1E6 (choose number in range of 10^5 to 10^8)
   * asymmetry: 0.001
   * m/z bin width: 1 (bin width of .1 versus 1 did not seem to make difference)
* **Mass detection (aka spectral filtering)
   * may want to use the centroid method since MZmine says data is already centroided
	* noise level = 100
   * assuming non-centroided data, wavelet transform looks good with following parameters:
	* noise level = 100
	* scale level = 10
	* window size (%) = 100
* **Chromatogram building (combine scans within each LCMS run):  
   Parameter values:
      * min time span (width) of peaks = .05 min (note: scans occur every .01 min)
      * min intensity of highest datapoint in peak = 500
      * mz tolerance (max diff btn mz values to be considered the same) = .01 (5ppm)
	* note: for a given scan, we have already choosen the centers of the peaks 
 	  to represent the mz value of that peak.  thus, we need to set mz tolerance
	  to the common spread of a peak, but rather the degree to which the center changes 
        *(I should verify that this is indeed true)
	* peaks appeared better looking when I used .01 as compared to higher values
* Chromatogram deconvolution: 
   * wavelets (XCMS) with the following parameters:
       * Signal/noise minimum threshold: 5
       * wavelet scales: .1 - 8 (a low min. value means narrow waves are not allowed)
       * peak duration range: .1 - 5 (I think this is similar to wavelet scales)
       * peak integration method: use smoothed data
* Deisotoping (parameter values inspired by a paper I came across): 
   * maximum charge of 6
   * retention time tolerance of 3 s
   * m/z tolerance of 0.03 Da
   * monotonic?  yes (I'm taking a guess on this, but I really don't know)
* Feature alignment across runs:
   * m/z tolerance of 0.05 Da
   * retention time tolerance of 1.5 min.
* Gap filling for missing value imputation: chose the "peak finding" method so that RT shifts can be accomodated
   * intensity tolerance: 20%
   * mz tolerance: .05 mz/ 5ppm
   * RT tolerance: .5 min
   * RT correction: yes
* Normalization
* Export peak list to csv for prediction analysis (to run in R)


## Option 2: run python code with pymzml package

* process_raw_data_and_do_prediction.py






