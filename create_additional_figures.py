

import os  #has several functions for manipulating files and directories
import sys
import pickle
import warnings
import time
import pdb #debugger
import math
import re #regular expression module

import numpy as np
import pandas as pd
print "Version of pandas: " , pd.__version__ #should be v0.16.1
from pandas.core.categorical import Categorical

import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt; plt.rcdefaults()
#import seaborn as sns #for prettier plots
#from matplotlib import style
#from matplotlib_style_utils import rstyle

import scipy as sp
from scipy import ndimage
from scipy import misc
from scipy import stats
from scipy.ndimage import imread
import scipy.ndimage as ndi

import sklearn
print "Version of sklearn: " , sklearn.__version__ #should be v0.16.1
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier, ExtraTreesClassifier
from sklearn import cross_validation
import sklearn.cross_validation as cv
from sklearn.cross_validation import cross_val_score, cross_val_predict
from sklearn import grid_search, metrics, datasets, preprocessing, feature_selection
from sklearn.feature_selection import SelectKBest
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn import svm, linear_model, neighbors, dummy, tree
from sklearn.pipeline import Pipeline
#from sklearn.linear_model import LogisticRegressionCV

sys.path.append('../SuPyLearner/supylearner')
import core #this is the main SuPyLearner code
from core import SuperLearner, cv_superlearner
from cross_val_utils import cross_val_predict_proba


np.random.seed(100)

preserveDir = "/home/ccotter/dengue_data_and_results_local/python_out/python_preserve/" #home PC
outDir = "/home/ccotter/dengue_data_and_results_local/python_out/" #home PC
#preserveDir = "/srv/scratch/ccotter/intermediate_data/python_preserve/" #mitra
#outDir = "/srv/scratch/ccotter/python_out/" 

def get_colors():
    """
        Return list of colors, defined by RGB values.
        Eventually may want to allow various parameters.
    """
    # These are the "Tableau 20" colors as RGB.  
        # They are ordered such that the first 10 will give you the tableau10 pallete
    tableau20 = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40),
                 (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127),
                 (188, 189, 34), (23, 190, 207),
                 (174, 199, 232), (255, 187, 120),  
                 (152, 223, 138),  (255, 152, 150),  
                  (197, 176, 213),  (196, 156, 148),  
                  (247, 182, 210),  (199, 199, 199),  
                  (219, 219, 141),  (158, 218, 229)] 
      
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
    for i in range(len(tableau20)):  
        r, g, b = tableau20[i]  
        tableau20[i] = (r / 255., g / 255., b / 255.)

    return tableau20

def main():

    ## Choose outcome variable ##
    outcome = "is.DEN"  
    #outcome = "is.DHF_DSS"

    ## Choose whether to exclude OFI patients ##
    NoOFI = False 

    ## Choose whether to exclude samples with initial DHF/DSS diagnosis ##
    NoInitialDHF = True 

    if outcome=="is.DEN":
        comparison_groups = "OFI vs. DENV using " #will appear in graph title
        FileNamePrefix = "OFI.v.DEN" #use all samples; exclude IR and PCR predictors
        NoInitialDHF = False #whether to exclude samples with initial DHF/DSS diagnosis
        NoOFI = False
    elif outcome=="is.DHF_DSS":
        if NoOFI == True:
            comparison_groups = "DF vs. DHF/DSS using " #will appear in graph title
            FileNamePrefix = "DF.v.DHFDSS"
        else:
            comparison_groups = "OFI/DF vs. DHF/DSS using " #will appear in graph title
            FileNamePrefix = "OFIDF.v.DHFDSS"

    ## Choose patient sample ##
    #patient_sample = "all"
    #patient_sample = "cohort_only"
    patient_sample = "hospital_only"

    figName = FileNamePrefix + '_combined_' + patient_sample + '.png'
    
    ## Combine results across various predictor sets
    predictor_desc_list = ["covarlist_all", "covarlist_noUltraX",
                            "covarlist_CohortRestrict","covarlist_genOnly"]

    ## Bar plot with bars ordered/grouped by algorithm and colors indicating predictors ##
    tableName = FileNamePrefix + '_' + "covarlist_all" + '_' + patient_sample + '.txt'
    resultsDF = pd.read_csv(preserveDir + 'R_' + tableName, sep=",")
    alg_names = resultsDF['Unnamed: 0'] #algorithm names
    initial_pos = np.arange(len(alg_names))*(len(predictor_desc_list)+.5)+len(predictor_desc_list)+.5
    bar_width = 1
    tableau20 = get_colors()
    #cycle through each predictor list
    plots = []
    for counter, predictor_desc in enumerate(predictor_desc_list):
        tableName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + '.txt'
        resultsDF = pd.read_csv(preserveDir + 'R_' + tableName, sep=",")
        measurements = np.array(resultsDF['cvAUC'])
        z = stats.norm.ppf(.95)
        SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
               ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
        alg_pos = initial_pos - counter 
        print "measurements: " , measurements
        print "alg_pos: " , alg_pos
        plot = plt.barh(bottom=alg_pos, width=measurements, height=bar_width,
                        xerr=SEs, align='center', alpha=1, 
                        color=tableau20[counter], label=predictor_desc)
        plots.append(plot)
    plt.xlabel = "cvAUC"
    plt.xlim(.5, 1)
    print "counter: " , counter
    plt.yticks(initial_pos - counter/2, alg_names)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outDir + figName)
    plt.show()            

if __name__ == '__main__':
    main()




