e###############################################################################
##################### Denge Dx Project Documentation ##########################
###############################################################################

Locations
=========

## Project code: Bitbucket repository
* https://bitbucket.org/cotterman/dengue-dx 
    * this one is to be retired, as bitbucket is improper location for non-code items
* https://bitbucket.org/cotterman/dengueDx-code
    * this one is active. As name implies, contains only code-related files
        * shortened name used with git commands: dcode
    * synced to the following folder:
        * on my PC: /home/carolyn/temp_Dengue_code (not maintained -- only for one-offs)
        * on Mitra: /users/ccotter/github_dengue


## Location of raw data (and other large dengue data files) on server: 
* wotan.lbl.gov: 
    * /srv/scratch/carolyn/dengue_project/raw_data (contains all raw .d files)
* amold.lbl.gov:
    * /srv/scratch/carolyn/mzml_serumR1, ../converted_serumR2 (contains serum mzML files)

* Note: The datasets I used for the R33 analysis are what Natalia put in folders labeled:
    * incoming/First batch 2012 (I relabled "/serumR1")
    * incoming/Second batch 2013 (I relabeled "/serumBadR2")
        * this data is bad and should probably not be used
    * incoming/Third batch 2014 (I relabeled "/serumR2)

## Location of other project materials:
* on my PC: /Dropbox/dengue_docs 
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

## Option 2: run python code with pymzml package

* process_raw_data_and_do_prediction.py