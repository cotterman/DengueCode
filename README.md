
Step 1: convert .d files to mzML format 
========================================

Note: msconvert works for many conversions, but for .d, one must use the Agilent DLLS, which only works on Windows.  


Step 2: Create features 
========================

## Option 1 - use MZmine2 software

* create batch queue file using MZmine GUI
* save this queue file as .xml 
* run MZmine from command line with .xml input option
* export XML file containing resulting peak list

## Option 2 - run python code with pymzml package

* process_raw_data.py

Step 3: Run prediction analysis 
================================

* run_prediction_in_python.py
