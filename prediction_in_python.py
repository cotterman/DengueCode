###############################################################################
###### Employ machine-learning algorithms in Python to predict Dengue  ########
###############################################################################


import os  #has several functions for manipulating files and directories
import pickle
import warnings
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy as sp
from scipy import ndimage
from scipy import misc
from scipy import stats
from scipy.ndimage import imread
import scipy.ndimage as ndi

from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
from sklearn.cross_validation import cross_val_score
from sklearn import grid_search
from sklearn import metrics
from sklearn import svm

from pandas.core.categorical import Categorical

inputsDir = "/home/carolyn/dengue_data_and_results_local/intermediate_data/"

VERBOSE = True
def log_statement(statement):
    if VERBOSE: print statement


def get_predictor_list(file_name, outcome):
    """
        Imports list of variables contained in file_name
        and replaces categorical variable names with binary equivalents
    """
    f = open(inputsDir + file_name, 'r')
    predictors = f.read().splitlines()

    #drop variables that are not appropriate for prediction type
    if outcome=="is.DEN":
        if predictors.count("PCR")>0 :
            predictors.remove("PCR")
        if predictors.count("IR")>0 :
            predictors.remove("IR")
    #drop variables that are not appropriate for prediction type
    if outcome=="is.DHF_DSS":
        

    #use binary equivalents of the categorical predictors
        #Torniquete -- less than 20 vs. 20+ = is.torniquete20plus
        #Pulso  -- strong or moderate vs. rapid or not palpable = is.pulse_danger
        #IR -- first vs. secondary infection (DENV positive only) = is.FirstInfection
        #PCR -- 1, 2, or 3 serotype (DENV positive only) = is.serotype1, is.serotype2, is.serotype3
    if predictors.count("Torniquete")>0 :
        predictors = [w.replace("Torniquete","is.torniquete20plus") for w in predictors]
    if predictors.count("Pulso")>0 :
        predictors = [w.replace("Pulso","is.pulse_danger") for w in predictors]
    if predictors.count("Sexo")>0 :
        predictors = [w.replace("Sexo","is.female") for w in predictors]
    if predictors.count("PCR")>0 :
        predictors.append("is.serotype1")
        predictors.append("is.serotype2")
        predictors.append("is.serotype3")
        predictors.remove("PCR")
    return(predictors)


def get_predictions_svm(data, outcome, predictors):
    #note: svm does not work with categorical data, only numeric and boolean data types

    clf = svm.SVC()
    clf.fit(data[predictors], data[outcome])
    fitted_values = clf.predict(data[predictors])
    return fitted_values

def get_performance_measures(actual, predicted):
    MSE = sum((actual - predicted)**2)
    return(MSE)

def convert_stings_to_categories(df):
    """
        Converts all variables of type "object" to type "category"
        of pandas dataframe
    """
    for col, dtype in zip(df.columns, df.dtypes):
        #print dtype, type(dtype)
        if dtype==np.object:
            print dtype            
            df[col]=df[col].astype("category")
    return df

def main():
    #filenames, outDir = parse_arguments() #filenames will be a list of (data, var_info)
    #os.chdir(outDir) #change pwd to output directory
    start_time_overall = time.time()

    #choose analysis type
    analysis = "Diagnose_OFI.v.DEN" #will use all samples; exclude IR and PCR predictors
    analysis = "Predict_OFIandDF.v.DHFandDSS" #exclude early DHF/DSS samples; exclude IR, PCR
    analysis = "Predict_DF.v.DHFandDSS" #exclude OFI and early DHF/DSS samples; use all predictors
    analysis = "Diagnose_DF.v.DHFandDSS" #all samples; most useful fancy predictors
    

    #obtain list of variables to use in prediction
    predictors = get_predictor_list("covarlist_all.txt", analysis)
    print "Predictors to use:\n" , predictors

    #create pandas dataframe with data that was cleaned in R
    df = pd.read_csv(inputsDir + "clin24_full_wImputedRF1.txt", sep='\t', 
        true_values=["True","Yes"], false_values=["False","No"], na_values=['NaN','NA']) 
    df[["is.female"]] = (df[["Sexo"]]=="female") #need this to be True/False binary
    df2 = df[predictors+[outcome]] #include only variables we will use
    #print "Number of rows in dataframe: " , df.shape[0], "\n"
    #print "Column names in dataframe: " , list(df2.columns.values), "\n"   

    #obtain cv-predictions using svm    
    fitted_values = get_predictions_svm(df2, outcome, predictors)
    print "Predicted values: " , fitted_values

    #obtain performance measures
    AUC = get_performance_measures(df2[outcome], fitted_values)
    print "cvAUC: " , AUC    

    log_statement("Total execution time, prediction program: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()
