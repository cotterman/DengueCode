###############################################################################
###### Employ machine-learning algorithms in Python to predict Dengue  ########
###############################################################################


import os  #has several functions for manipulating files and directories
import sys
sys.path.append('../SuPyLearner/supylearner')
import core #this is the main SuPyLearner code
from core import SuperLearner, cv_superlearner

import pickle
import warnings
import time

import numpy as np
import pandas as pd
from pandas.core.categorical import Categorical
import matplotlib.pyplot as plt

import scipy as sp
from scipy import ndimage
from scipy import misc
from scipy import stats
from scipy.ndimage import imread
import scipy.ndimage as ndi

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn import cross_validation
from sklearn.cross_validation import cross_val_score
from sklearn import grid_search, metrics, datasets
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn import svm, linear_model, neighbors
#from sklearn.linear_model import LogisticRegressionCV

np.random.seed(100)

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
    #if outcome=="is.DHF_DSS":
        

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

    clf = svm.SVC(probability=True)
    clf.fit(data[predictors], data[outcome])
    fitted_values = clf.predict(data[predictors])
    return fitted_values

def get_performance_measures(actual, predicted):
    MSE = sum((actual.astype(int) - predicted)**2)
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

    #choose analysis type (to-do)
    analysis = "Diagnose_OFI.v.DEN" #will use all samples; exclude IR and PCR predictors
    analysis = "Predict_OFIandDF.v.DHFandDSS" #exclude early DHF/DSS samples; exclude IR, PCR
    analysis = "Predict_DF.v.DHFandDSS" #exclude OFI and early DHF/DSS samples; use all predictors
    analysis = "Diagnose_DF.v.DHFandDSS" #all samples; most useful fancy predictors
    outcome = "is.DEN"

    #obtain list of variables to use in prediction
    predictors = get_predictor_list("covarlist_all.txt", outcome)
    print "Predictors to use:\n" , predictors

    #create pandas dataframe with data that was cleaned in R
    df = pd.read_csv(inputsDir + "clin24_full_wImputedRF1.txt", sep='\t', 
        true_values=["True","Yes"], false_values=["False","No"], na_values=['NaN','NA']) 
    #df["is.DEN"] = df["is.DEN"].astype(int) #ML functions give predicted probs only for int outcomes
    df[["is.female"]] = (df[["Sexo"]]=="female") #need this to be True/False binary
    print "Number of rows in dataframe: " , df.shape[0], "\n"
    #print "Column names in dataframe: " , list(df.columns.values), "\n"  
    X = df[predictors] #include only variables we will use 
    y = df[outcome]
    #X2, y2=datasets.make_classification(n_samples=88, n_features=95) #generate toy data for testing

    #run Super Learner
    RF = RandomForestClassifier()
    #logitL1=linear_model.LogisticRegressionCV(penalty='l1', solver='liblinear') 
        #todo: explore l1 path options for optimizing penalty coefficient (e.g., l1_min_c)
    #logitL2=linear_model.RidgeClassifierCV(normalize=True) 
    #nn=neighbors.KNeighborsClassifier(n_neighbors=4) #Classes are ordered by lexicographic order.
    #gradB_dev=GradientBoostingClassifier(loss='deviance') #for probabilistic outputs
    adaBoost=AdaBoostClassifier() 
    #myLDA = LDA(solver='lsqr', shrinkage='auto') #optimal shrinkage is calculated analytically
    #myQDA = QDA() 
    #svm1=svm.SVC(kernel='rbf', probability=True)     
    #svm2=svm.SVC(kernel='poly', probability=True)
    #svmL2=svm.LinearSVC(penalty='l2') #penalized linear support vector classification
    #todo: spectral clustering classifier?
    #todo: write a means classifier myself
    lib=[svm1, svm2]
    libnames=["SVM-rbf", "SVM-poly"]
    #lib=[RF, logitL1, logitL2, nn, gradB_dev, adaBoost, myLDA, myQDA, svm1, svm2, svmL2]
    #libnames=["RF", "LogitL1", "LogitL2", "Nearest Neighbor","Gradient Boosting",
    #           "AdaBoost", "LDA", "QDA", "SVM-rbf", "SVM-poly", "SVM-L2"]

    sl=SuperLearner(lib, libnames, loss="nloglik", K=2, stratifyCV=True)
    sl.fit(X, y)
    sl.summarize()
    cv_superlearner(sl, X, y, K=2)

    log_statement("Total execution time, prediction program: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()
