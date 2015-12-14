

import os  #has several functions for manipulating files and directories
import sys
import pickle
import warnings
import time
import pdb #debugger
import math
import re #regular expression module
from collections import namedtuple

import numpy as np
import pandas as pd
print "Version of pandas: " , pd.__version__ #should be v0.16.1
from pandas.core.categorical import Categorical

import matplotlib as mpl
mpl.use('Agg') #avoids running an X server (to avoid errors with remote runs)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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


sys.path.append('../SuPyLearner/supylearner')
import core #this is the main SuPyLearner code
from core import SuperLearner, cv_superlearner
from cross_val_utils import cross_val_predict_proba


np.random.seed(100)

#outDir = "/home/ccotter/dengue_data_and_results_local/python_out/" #home PC
outDir = "/users/ccotter/python_out/" #nandi

GraphInfo = namedtuple("GraphInfo", 
    ["patient_sample","comparison_groups", "FileNamePrefix", "NoInitialDHF", "NoOFI"])


def get_colors():
    """Return list of colors, defined by RGB values.
       
    """
    #Eventually may want to allow various parameters.
    
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


def parse_args(patient_sample, outcome, NoOFI, NoInitialDHF):
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
    return GraphInfo(patient_sample, comparison_groups, FileNamePrefix, NoInitialDHF, NoOFI)

def create_LCMS_topJ_barplot(ginfo):
  
    ## Set graph features
    fig, ax = plt.subplots()
    bar_width = 1
    #palette will be ordered from light to dark by default
    colors_noLCMS = sns.light_palette("charcoal", input="xkcd", reverse=True, n_colors=4) 
    colors_ttest = sns.light_palette("forest green", input="xkcd", reverse=True, n_colors=4) 
    colors_topRF = sns.light_palette("dark blue", input="xkcd", reverse=True, n_colors=4) 
    plt.figure(figsize=(6.7,8)) 

    ## Decide which data to include
    temp = range(7) + [744]
    LCMSvars_list = temp[::-1] #want 0 to be high on graph - plot it last
    subsetAlg_list = ['ttest','topRF']
    method_list = ['M3','M2','M1'] #M1 will be high on graph (good)

    ## Prepare data to be graphed    
    ttestDF = pd.read_csv(outDir + 'R_OFIDF.v.DHFDSS_BestSubset_ttest_NPserum.txt', sep=",")
    ttestDF['mlabel'] = ttestDF['Unnamed: 0']
    topRFDF = pd.read_csv(outDir + 'R_OFIDF.v.DHFDSS_BestSubset_topRF_NPserum.txt', sep=",")
    topRFDF['mlabel'] = topRFDF['Unnamed: 0']

    ## To fill in during loop
    alg_pos = []
    measurements = []
    colors = []
    method_labels = []
    colors_legend = []
    ytick_labels = []
    ytick_positions = []
    ymax = 0

    ## Cycle through for number of LCMS vars
    for counter, LCMSvars in enumerate(LCMSvars_list):
        ytick_labels.append("Clinical + \n"+str(LCMSvars)+" LCMS features")
        if LCMSvars==0:
            ytick_positions.append(ymax+1)
        else:
            ytick_positions.append(ymax+3)

        #cycle through LCMS feature selection method ('ttest' and 'topRF')
        for num, subsetAlg in enumerate(subsetAlg_list):

            #cycle through methods for treating clinical info (1, 2, and 3)
            for method in method_list:
                #only do once when not adding any LCMS features
                    #does not make sense to compare the ttest and topRF subsetAlg approaches
                if LCMSvars == 0 and num>0:
                    continue
                rowlab = method + '_LCMS_' + str(LCMSvars)
                if subsetAlg=="ttest": 
                    df = ttestDF
                    colors_list = colors_ttest
                elif subsetAlg=="topRF": 
                    df = topRFDF
                    colors_list = colors_topRF
                if LCMSvars == 0:
                    colors_list = colors_noLCMS
                    method_labels.append(method)
                    colors_legend.append(colors_list[int(method[1])-1])
                elif LCMSvars == 1:
                    method_labels.append(subsetAlg + ' ' + method)
                    colors_legend.append(colors_list[int(method[1])-1])
                measurements.append(float(df.loc[df['mlabel']==rowlab]['cvAUC']))
                alg_pos.append(ymax)
                #color list was determined by subsetAlg; saturation is based on method
                colors.append(colors_list[int(method[1])-1])
                ymax += bar_width
        #add space between groups of bars segmented by number of LCMS features
        ymax += bar_width

    plt.barh(bottom=alg_pos, width=measurements, height=bar_width,
                    align='center', alpha=1, color=colors)
    plt.yticks(ytick_positions, ytick_labels) #size=16
    plt.xlim(.5, 1)
    plt.ylim(-2, ymax+6) 
    #make left spacing large enough for labels.  Default is  .1, .9, .9, .1
    plt.subplots_adjust(left=.22, right=.9, top=.9, bottom=.1)
    lhandles = []
    for mycolor in colors_legend[::-1]:
        hand = mpatches.Patch(color=mycolor)
        lhandles.append(hand)
    leg = plt.legend((lhandles), (method_labels[::-1]), ncol=3)
    #plt.tight_layout()

    plt.savefig(outDir + 'LCMS_subsets_OFIDF.v.DHFDSS.eps', dpi=1200) 
    plt.close()  

def create_LCMS_barplot(ginfo, LCMScompare):
    """Bar plot with bars grouped by predictor set and colors indicating LCMS run  

        LCMScompare = "NPbins_v_RPbins" to compare NP vs. RP using binned data
        LCMScompare = "NPbins_v_MassHuntNP" to comapre NP binned vs. NP mass hunter
    """

    if LCMScompare == "NPbins_v_RPbins":
        # LCMS data to loop through
        inLCMSData_list = ['NPbins50x50', 'RPbins50x50']
        # labels to appear in graph legend
        inLCMSData_desc = ['Normal phase, 50x50 intensity grid',
                       'Reverse phase, 50x50 intensity grid'] 
        #xkcd colors
        colors = ["light burgundy", "medium blue"] 

    elif LCMScompare == "NPbins_v_MassHuntNP":
        # LCMS data to loop through
        inLCMSData_list = ['NPbins50x50', 'MassHuntNP']
        # labels to appear in graph legend
        inLCMSData_desc = ['Normal phase, 50x50 intensity grid',
                       'Normal phase, Mass Hunter'] 
        #xkcd colors
        colors = ["light burgundy", "tangerine"] 

    elif LCMScompare == "NonInvasives":
        # LCMS data to loop through
        inLCMSData_list = ['SalivaMH','UrineMH']
        # labels to appear in graph legend
        inLCMSData_desc = ['Saliva','Urine'] 
        #xkcd colors
        colors = ["dark lilac","teal green"]

    figName = ginfo.FileNamePrefix + '_' + LCMScompare
    #first name listed will appear closest to bottom of y-axis
    predcat_names = ['Clinical+LCMS','LCMS only','Clinical only'] 
    predictor_desc = "covarlist_all"
    #will appear in alphabetical order
    alg_list = ['Gradient Boost','AdaBoost','Random Forests', 'Super Learner'] #
    initial_pos = np.arange(len(predcat_names)*len(alg_list))*(
        len(inLCMSData_list)+1)+len(inLCMSData_list)+1
    #initial_pos = np.arange(len(predcat_names)*len(alg_list)+len(predcat_names))*(
    #    len(inLCMSData_list)+1)+len(inLCMSData_list)+1
    bar_width = 1
    mycolors = sns.xkcd_palette(colors)
    plt.figure(figsize=(6.7,8))
    #cycle through each inLCMSData value
    plots = []
    for counter, inLCMSData in enumerate(inLCMSData_list):
        resultsDF = pd.DataFrame()
        for myc, predcat in enumerate(predcat_names):            
            if predcat=='Clinical only':
                tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_'+inLCMSData+'patients'
            elif predcat=='LCMS only':
                tableName = ginfo.FileNamePrefix + '_covarlist_' + inLCMSData
            elif predcat=='Clinical+LCMS':
                tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_' + inLCMSData
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
        print "tablename: " , tableName 
        print "resultsDF: " , resultsDF
        print "measurements: " , measurements
        print "ypos: " , ypositions
        plot = plt.barh(bottom=ypositions, width=measurements, height=bar_width,
                        xerr=SEs, error_kw=dict(ecolor='.1', lw=1, capsize=1, capthick=1),
                        align='center', alpha=1, 
                        color=mycolors[counter], label=inLCMSData_desc[counter])
        #add numeric values to plot
        xpos = np.array(resultsDF['ci_low']) -.05
        ypos = ypositions - .3
        mytext = ["%.2f" % x for x in measurements]    
        for place, text in enumerate(mytext):
            plt.text(xpos[place], ypos[place], text, color="white", fontsize=12)
        plots.append(plot)

    #add labels for algorithms
    plt.yticks(initial_pos - counter/2 - .5, alg_list*len(predcat_names))

    #add labels for predictors used (LCMS vs clinical vs combo)
    print "len(predcat_names): " , len(predcat_names)
    stepsize = len(predcat_names)*len(alg_list)
    for myc, pred_name in enumerate(predcat_names):
        ypos = stepsize*(myc+1) + .6
        plt.text(.5, ypos, pred_name, color="black", fontsize=14)

    plt.xlabel = "cvAUC"
    plt.xlim(.5, 1)
    plt.ylim(0,max(initial_pos)+4)
    print "counter: " , counter
    plt.legend(ncol=1)
    plt.tight_layout()
    plt.savefig(outDir + figName + '.eps', dpi=1200)
    plt.close() 


def create_predSets_barplot(ginfo):
  
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
    figName = ginfo.FileNamePrefix + '_combined_' + ginfo.patient_sample + '_PC'
    tableName = ginfo.FileNamePrefix + '_' + "covarlist_all" + '_' + ginfo.patient_sample + '.txt'
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
        tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_' + ginfo.patient_sample + '.txt'
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

def create_impDum_barplot(ginfo):    
    #runs to loop through
    suffix_list = ["_impDums", "", "_dumsOnly"]
    #labels to appear in graph legend
    list_desc = ["Clinical values + imputation indicators",
                "Clinical values only",
                "Imputation indicators only"]
    predictor_desc = "covarlist_all" 
    figName   = ginfo.FileNamePrefix + '_' + predictor_desc + ginfo.patient_sample + '_ImpAnalysis' 
    tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_' + ginfo.patient_sample + '.txt'
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
        tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_' + ginfo.patient_sample + suffix + '.txt'
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

def create_testData_barplot(ginfo):
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
    figName   = ginfo.FileNamePrefix + '_' + predictor_desc + '_TestSets' 
    tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_cohortTest' + '.txt'
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
    for counter, ginfo.patient_sample in enumerate(patient_list):
        tableName = ginfo.FileNamePrefix + '_' + predictor_desc + '_' + ginfo.patient_sample + '.txt'
        resultsDF = pd.read_csv(outDir + 'R_' + tableName, sep=",")
        #results from test set do not come with CIs
        if ginfo.patient_sample not in ['hospitalTest','cohortTest']:
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


def main():

    ## Choose outcome variable 
    outcome = "is.DEN"  
    #outcome = "is.DHF_DSS"

    ## Choose whether to exclude OFI patients 
    NoOFI = False #only applicable for is.DHF_DSS

    ## Choose whether to exclude samples with initial DHF/DSS diagnosis 
    NoInitialDHF = True #only applicable for is.DHF_DSS

    ## Choose patient sample 
    patient_sample = "all"
    #patient_sample = "cohort_only"
    #patient_sample = "hospital_only"

    ## Obtain graph title, filename prefix, etc.
    ginfo = parse_args(patient_sample, outcome, NoOFI, NoInitialDHF)

    ## Bart plot which compares methods to choose best subset of LCMS features
    create_LCMS_topJ_barplot(ginfo)

    ## Bar plot with bars grouped by predictor set and colors indicating LCMS run 
    #LCMScompare = "NonInvasives" #"NonInvasives", "NPbins_v_MassHuntNP", "NPbins_v_RPbins"
    #create_LCMS_barplot(ginfo, LCMScompare)   

    ## Bar plot with bars ordered/grouped by algorithm and colors indicating predictors sets 
    #create_predSets_barplot(ginfo)

    ## Bar plot with bars grouped by algorithm and colors indicating imputation dummy inclusion 
    #create_impDum_barplot(ginfo)  

    ## Bar plot with bars grouped by algorithm and colors indicating patient sample 
    #create_testData_barplot(ginfo)

if __name__ == '__main__':
    main()




