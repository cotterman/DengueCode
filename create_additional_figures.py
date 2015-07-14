

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

import matplotlib as mpl
mpl.use('Agg') #avoids running an X server (to avoid errors with remote runs)
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt; plt.rcdefaults()
#import seaborn as sns #for prettier plots
#from matplotlib import style
#from matplotlib_style_utils import rstyle
import matplotlib.lines as mlines
import seaborn as sns

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

#inputsDir = "/home/ccotter/dengue_data_and_results_local/intermediate_data/" #home PC
#outDir = "/home/ccotter/dengue_data_and_results_local/python_out/" #home PC
inputsDir = "/srv/scratch/ccotter/intermediate_data/" #mitra and nandi
outDir = "/users/ccotter/python_out/" #mitra and nandi

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

    ## Choose which plots to create ##
    create_master_barplot = True
    create_testData_barplot = True 

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

    ## Bar plot with bars ordered/grouped by algorithm and colors indicating predictors ##
    if create_master_barplot==True:
        #combine results across these predictor sets
        predictor_desc_list = ["covarlist_all", "covarlist_noUltraX",
                            "covarlist_CohortRestrict","covarlist_genOnly"]
        pred_label_list = ['Basic+Lab+Costly','Basic+Lab','Clinic Collected','Basic Only']
        figName = FileNamePrefix + '_combined_' + patient_sample
        tableName = FileNamePrefix + '_' + "covarlist_all" + '_' + patient_sample + '.txt'
        resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
        alg_names = resultsDF['Unnamed: 0'] #algorithm names
        #could consider linspace instead of np.arange
        initial_pos = np.arange(len(alg_names))*(
            len(predictor_desc_list)+.5)+len(predictor_desc_list)+.5
        bar_width = 1
        mycolors = sns.cubehelix_palette(4, reverse=True) #same hue with diff intensities
        plt.figure(figsize=(6.7,8))
        #cycle through each predictor list
        plots = []
        for counter, predictor_desc in enumerate(predictor_desc_list):
            tableName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + '.txt'
            resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
            measurements = np.array(resultsDF['cvAUC'])
            z = stats.norm.ppf(.95)
            SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
                   ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
            alg_pos = initial_pos - counter 
            print "measurements: " , measurements
            print "alg_pos: " , alg_pos
            plot = plt.barh(bottom=alg_pos, width=measurements, height=bar_width,
                            xerr=SEs, error_kw=dict(ecolor='.1', lw=1, capsize=1.5, capthick=1),
                            align='center', alpha=1, 
                            color=mycolors[counter], label=pred_label_list[counter])
            plots.append(plot)
        plt.xlabel = "cvAUC"
        plt.xlim(.5, 1)
        plt.ylim(0,max(initial_pos)+1) 
        print "counter: " , counter
        plt.yticks(initial_pos - counter/2, alg_names)
        plt.legend(title="Clinical Variables")
        plt.tight_layout()
        plt.savefig(outDir + figName + '.eps', dpi=1200)
        plt.close()           

    ## Bar plot with bars grouped by algorithm and colors indicating patient sample ##
    if create_testData_barplot==True:
        # patient samples to loop through
        patient_list = ['all_studyDum', 'all_noDums',
                        'hospital_only','hospitalTest',
                        'cohort_only', 'cohortTest'] 
                        #,, 
        # labels to appear in graph legend
        patient_desc = ['Hospital + clinic with study indicator (CV results)',
                        'Hospital + clinic without study indicator (CV results)',
                        'Hospital only (CV results)',
                        'Predictions for hospital after fitting to clinic',
                        'Clinic only (CV results)',
                        'Predictions for clinic after fitting to hospital']
        predictor_desc = "covarlist_CohortRestrict"
        figName   = FileNamePrefix + '_' + predictor_desc + '_TestSets' 
        tableName = FileNamePrefix + '_' + predictor_desc + '_cohortTest' + '.txt'
        resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
        alg_names = resultsDF['Unnamed: 0'] #algorithm names
        print "alg_names: " , alg_names
        initial_pos = np.arange(len(alg_names))*(
            len(patient_list)+1)+len(patient_list)+1
        bar_width = 1
        #hues inspired by list(reversed(sns.color_palette("Paired",12)))
        mycolors = [(0.42485198495434734, 0.2511495584950722, 0.60386007743723258),
                    (0.78329874347238004, 0.68724338552531095, 0.8336793640080622),
                    (0.12572087695201239, 0.47323337360924367, 0.707327968232772),
                    (0.65098041296005249, 0.80784314870834351, 0.89019608497619629),
                    (0.89059593116535862, 0.10449827132271793, 0.11108035462744099),
                    (0.98320646005518297, 0.5980161709820524, 0.59423301088459368)]
        plt.figure(figsize=(6.7,8))
        #cycle through each patient list
        plots = []
        for counter, patient_sample in enumerate(patient_list):
            tableName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + '.txt'
            resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
            #results from test set do not come with CIs
            if patient_sample not in ['hospitalTest','cohortTest']:
                measurements = np.array(resultsDF['cvAUC'])
                z = stats.norm.ppf(.95)
                SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
                       ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
            else:
                measurements = np.array(resultsDF['AUC'])
            alg_pos = initial_pos - counter 
            print "measurements: " , measurements
            print "alg_pos: " , alg_pos
            plot = plt.barh(bottom=alg_pos, width=measurements, height=bar_width,
                            xerr=SEs, error_kw=dict(ecolor='.1', lw=1, capsize=1, capthick=1),
                            align='center', alpha=1, 
                            color=mycolors[counter], label=patient_desc[counter])
            plots.append(plot)
        plt.xlabel = "cvAUC"
        plt.xlim(.5, 1)
        plt.ylim(0,max(initial_pos)+2)
        print "counter: " , counter
        plt.yticks(initial_pos - counter/2, alg_names)
        plt.legend(prop={'size':6})
        plt.tight_layout()
        plt.savefig(outDir + figName + '.eps', dpi=1200)
        plt.close()            

if __name__ == '__main__':
    main()




