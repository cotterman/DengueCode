

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

inputsDir = "/home/ccotter/dengue_data_and_results_local/intermediate_data/" #home PC
outDir = "/home/ccotter/dengue_data_and_results_local/python_out/" #home PC
#inputsDir = "/srv/scratch/ccotter/intermediate_data/" #mitra and nandi
#outDir = "/users/ccotter/python_out/" #mitra and nandi

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
    create_predSets_barplot = False
    create_impDum_barplot = False
    create_testData_barplot = False 
    create_LCMS_barplot = True

    ## Choose outcome variable ##
    outcome = "is.DEN"  
    #outcome = "is.DHF_DSS"

    ## Choose whether to exclude OFI patients ##
    NoOFI = True #only applicable for is.DHF_DSS

    ## Choose whether to exclude samples with initial DHF/DSS diagnosis ##
    NoInitialDHF = True #only applicable for is.DHF_DSS

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

    ## Bar plot with bars grouped by predictor set and colors indicating LCMS run (RP/NP) ##
    if create_LCMS_barplot==True:
        # LCMS data to loop through
        inLCMSData_list = ['NPbins50x50', 'RPbins50x50']
        # labels to appear in graph legend
        inLCMSData_desc = ['Normal phase, 50x50 intensity grid',
                           'Reverse phase, 50x50 intensity grid']  
        figName = FileNamePrefix + '_LCMScompare'
        #first name listed will appear closest to bottom of y-axis
        predcat_names = ['Clinical+LCMS','LCMS only','Clinical only'] 
        predictor_desc = "covarlist_all"
        #will appear in alphabetical order
        alg_list = ['Gradient Boost','AdaBoost','Random Forests','Super Learner']
        initial_pos = np.arange(len(predcat_names)*len(alg_list))*(
            len(inLCMSData_list)+1)+len(inLCMSData_list)+1
        #initial_pos = np.arange(len(predcat_names)*len(alg_list)+len(predcat_names))*(
        #    len(inLCMSData_list)+1)+len(inLCMSData_list)+1
        bar_width = 1
        colors = ["medium blue","light burgundy"]
        mycolors = sns.xkcd_palette(colors)
        #plt.figure(figsize=(4,6))
        #cycle through each inLCMSData value
        plots = []
        for counter, inLCMSData in enumerate(inLCMSData_list):
            resultsDF = pd.DataFrame()
            for myc, predcat in enumerate(predcat_names):            
                if predcat=='Clinical only':
                    tableName = FileNamePrefix + '_' + predictor_desc + '_'+inLCMSData+'patients'
                elif predcat=='LCMS only':
                    tableName = FileNamePrefix + '_covarlist_' + inLCMSData
                elif predcat=='Clinical+LCMS':
                    tableName = FileNamePrefix + '_' + predictor_desc + '_' + inLCMSData
                rDF = pd.read_csv(outDir + 'R_' + tableName + '.txt', sep=",")
                #sort by alg_list so numbers match with labels
                rDF.sort(columns=['Unnamed: 0'], axis=0, ascending=False, inplace=True)
                prelim = rDF[rDF['Unnamed: 0'].isin(alg_list)]
                resultsDF = pd.concat([resultsDF,prelim],axis=0)
            alg_list.sort(reverse=True) #to match with numbers
            measurements = np.array(resultsDF['cvAUC'])
            xpositions = np.array(resultsDF['ci_low'])
            z = stats.norm.ppf(.95)
            SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
                   ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
            ypositions = initial_pos - counter  
            print "resultsDF: " , resultsDF
            print "measurements: " , measurements
            print "ypos: " , ypositions
            plot = plt.barh(bottom=ypositions, width=measurements, height=bar_width,
                            xerr=SEs, error_kw=dict(ecolor='.1', lw=1, capsize=1, capthick=1),
                            align='center', alpha=1, 
                            color=mycolors[counter], label=inLCMSData_desc[counter])
            #add numeric values to plot
            xpos = np.array(resultsDF['ci_low']) -.02
            ypos = ypositions - .3
            mytext = ["%.2f" % x for x in measurements]    
            for place, text in enumerate(mytext):
                plt.text(xpos[place], ypos[place], text, color="white")
            plots.append(plot)

        #add labels for algorithms
        plt.yticks(initial_pos - counter/2 - .5, alg_list*len(predcat_names))

        #add labels for predictors used (LCMS vs clinical vs combo)
        print "len(predcat_names): " , len(predcat_names)
        stepsize = len(predcat_names)*len(alg_list)
        for myc, pred_name in enumerate(predcat_names):
            ypos = stepsize*(myc+1) + .6
            plt.text(.5, ypos, pred_name, color="black")

        plt.xlabel = "cvAUC"
        plt.xlim(.5, 1)
        plt.ylim(0,max(initial_pos)+4)
        print "counter: " , counter
        plt.legend()
        plt.tight_layout()
        plt.savefig(outDir + figName + '.png', dpi=1200)
        plt.close()    


    ## Bar plot with bars ordered/grouped by algorithm and colors indicating predictors sets ##
    if create_predSets_barplot==True:
      
        sns.set_style("whitegrid") #necessary for getting back graph frame
        mpl.rcParams['lines.color'] = 'white'
        mpl.rcParams['text.color'] = 'white'
        fig, ax = plt.subplots()

        #combine results across these predictor sets
        #predictor_desc_list = ["covarlist_all", "covarlist_noUltraX",
        #                    "covarlist_CohortRestrict","covarlist_genOnly"]
        predictor_desc_list = ["covarlist_all", "covarlist_noUltraX","covarlist_genOnly"] #JSM
        #pred_label_list = ['Basic+Lab+Costly','Basic+Lab','Clinic Collected','Basic Only']
        pred_label_list = ['Basic+Lab+Costly','Basic+Lab','Basic Only'] #JSM
        figName = FileNamePrefix + '_combined_' + patient_sample + '_PC'
        tableName = FileNamePrefix + '_' + "covarlist_all" + '_' + patient_sample + '.txt'
        resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
        alg_names = resultsDF['Unnamed: 0'] #algorithm names
        #alg_names = ['Super Learner'] #for shortened graph (JSM poster)
        #could consider linspace instead of np.arange
        initial_pos = np.arange(len(alg_names))*(
            len(predictor_desc_list)+.5)+len(predictor_desc_list)+.5
        bar_width = 1
        #mycolors = sns.cubehelix_palette(4, reverse=True, dark=.3, rot=-.4) #rot=-.4 for green
        mycolors = sns.cubehelix_palette(4, start=3, dark=.3, reverse=True) #JSM poster (pink)
        plt.figure(figsize=(8,8)) #dissertation: (6.7,8), poster: (6,7) and (6,2)
        plt.grid(b=True, which='both', axis='both',color='0.3',linestyle='-')
        
        #cycle through each predictor list
        plots = []
        for counter, predictor_desc in enumerate(predictor_desc_list):
            tableName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + '.txt'
            print "tableName: " , tableName
            resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
            #rDF = pd.read_csv(outDir + 'R_' + tableName, sep=",") #JSM poster
            #resultsDF = rDF[rDF['Unnamed: 0'].isin(['Super Learner'])] #JSM poster
            measurements = np.array(resultsDF['cvAUC'])
            z = stats.norm.ppf(.95)
            SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
                   ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
            alg_pos = initial_pos - counter 
            print "measurements: " , measurements
            print "alg_pos: " , alg_pos
            plot = plt.barh(bottom=alg_pos, width=measurements, height=bar_width,
                            xerr=SEs, error_kw=dict(ecolor='1', lw=1, capsize=1.5, capthick=1),
                            align='center', alpha=1, 
                            color=mycolors[counter], label=pred_label_list[counter])
            plots.append(plot)
        plt.xlabel("cvAUC",color='white', size=16)
        plt.xlim(.5, 1)
        plt.ylim(0,max(initial_pos)+1) 
        print "counter: " , counter
        plt.yticks(initial_pos - counter/2, alg_names, color='white', size=16)
        #plt.yticks(np.arange(len(predictor_desc_list))+1, pred_label_list, color='white') #JSM poster
        plt.xticks(color='white', size=16)
        leg = plt.legend(title="Clinical Variables")
        plt.legend(prop={'size':16})
        for text in leg.get_texts():
            plt.setp(text, color='white', size=16)
        #plt.title('Performance for distinguishing OFI from dengue')
        plt.tight_layout()
        plt.savefig(outDir + figName + '_talk.png', dpi=1200, transparent=True) #transparent=True for ppt
        plt.close()     

    ## Bar plot with bars grouped by algorithm and colors indicating imputation dummy inclusion ##  
    if create_impDum_barplot==True:    
        #runs to loop through
        suffix_list = ["_impDums", "", "_dumsOnly"]
        #labels to appear in graph legend
        list_desc = ["Clinical values + imputation indicators",
                    "Clinical values only",
                    "Imputation indicators only"]
        predictor_desc = "covarlist_all" 
        figName   = FileNamePrefix + '_' + predictor_desc + patient_sample + '_ImpAnalysis' 
        tableName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + '.txt'
        resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
        alg_names = resultsDF['Unnamed: 0'] #algorithm names
        print "alg_names: " , alg_names
        initial_pos = np.arange(len(alg_names))*(
            len(suffix_list)+1)+len(suffix_list)+1
        bar_width = 1
        colors = ["amber","windows blue","greyish"]
        mycolors = sns.xkcd_palette(colors)
        plt.figure(figsize=(6.7,8))
        #cycle through each patient list
        plots = []
        for counter, suffix in enumerate(suffix_list):
            tableName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + suffix + '.txt'
            resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
            measurements = np.array(resultsDF['cvAUC'])
            z = stats.norm.ppf(.95)
            SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
                   ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
            alg_pos = initial_pos - counter 
            print "measurements: " , measurements
            print "alg_pos: " , alg_pos
            plot = plt.barh(bottom=alg_pos, width=measurements, height=bar_width,
                            xerr=SEs, error_kw=dict(ecolor='.1', lw=1, capsize=1, capthick=1),
                            align='center', alpha=1, 
                            color=mycolors[counter], label=list_desc[counter])
            plots.append(plot)
        plt.xlabel = "cvAUC"
        plt.xlim(.5, 1)
        plt.ylim(0,max(initial_pos)+2)
        print "counter: " , counter
        plt.yticks(initial_pos - counter/2, alg_names)
        plt.legend(prop={'size':8})
        plt.tight_layout()
        plt.savefig(outDir + figName + '.eps', dpi=1200)
        plt.close()    


    ## Bar plot with bars grouped by algorithm and colors indicating patient sample ##
    if create_testData_barplot==True:
        # patient samples to loop through
        patient_list = ['all_studyDum', 'all',
                        'hospital_only','hospitalTest',
                        'cohort_only', 'cohortTest'] 
                        #,, 
        # labels to appear in graph legend
        patient_desc = ['Hospital + clinic with study indicator (CV results)',
                        'Hospital + clinic without study indicator (CV results)',
                        'Hospital only (CV results)',
                        'Predictions for hospital after training with clinic',
                        'Clinic only (CV results)',
                        'Predictions for clinic after training with hospital']
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




