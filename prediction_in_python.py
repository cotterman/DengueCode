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
print "Version of pandas: " , pd.__version__ #should be v0.16.1
from pandas.core.categorical import Categorical
import matplotlib.pyplot as plt

import scipy as sp
from scipy import ndimage
from scipy import misc
from scipy import stats
from scipy.ndimage import imread
import scipy.ndimage as ndi

import sklearn
print "Version of sklearn: " , sklearn.__version__ #should be v0.16.1
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn import cross_validation
from sklearn.cross_validation import cross_val_score, cross_val_predict
from sklearn import grid_search, metrics, datasets
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn import svm, linear_model, neighbors
#from sklearn.linear_model import LogisticRegressionCV

import pdb #debugger

np.random.seed(100)

#inputsDir = "/home/carolyn/dengue_data_and_results_local/intermediate_data/" #home PC
#inputsDir = "/home/nboley/Desktop/dengue_data_and_results_local/intermediate_data/" #home PC, diff login
inputsDir = "/srv/scratch/ccotter/intermediate_data/" #mitra

class RidgeClassifier(linear_model.RidgeClassifier):
    def predict_proba(self, X): 
        #pdb.set_trace()
        coefs = self.coef_
        print type(coefs), coefs
        z = np.dot(X, np.transpose(coefs))
        print type(z)
        print z.shape
        # return the predicted probabilities (of getting outcome 1)
        predict_proba = 1.0 / (1.0 + np.exp(-z))
        return predict_proba[:,0]

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


def get_predictions_svm(X, y):
    #note: svm does not work with categorical data, only numeric and boolean data types

    clf = svm.SVC(probability=True)
    clf.fit(X, y)
    fitted_values = clf.predict_proba(X)[:, 0]
    #preds_cv = cross_val_predict(clf, X, y)
    print clf.classes_ #returns array containing order classes appear in fitted_values (todo)
    return fitted_values

def measure_performance(y, pred_prob, threshold):
    """
        Obtain measures that depend on specified classification threshold
    """
    pred_y = (pred_prob>threshold).astype(int) #turn predicted probs into 0/1 predicted value
    #counts to use in performance measures
    true_positives_count = (np.absolute(y + pred_y - 2) < 1e-6).astype(int).sum()
    positives_count = (np.absolute(y-1)<1e-6).astype(int).sum()
    predicted_positives_count = sum(pred_y==1)
    true_negatives_count = ((np.absolute(y + pred_y)) < 1e-6).astype(int).sum()
    negatives_count = (np.absolute(y)<1e-6).astype(int).sum()
    predicted_negatives_count = sum(pred_y==0)
    error_count = (np.absolute(pred_y - y) > 1e-6).astype(int).sum()
    #calculate performance measures
    sensitivity = true_positives_count.astype(float)/ positives_count
    specificity = true_negatives_count.astype(float) / negatives_count
    PPV = true_positives_count.astype(float) / predicted_positives_count 
    NPV = true_negatives_count.astype(float) / predicted_negatives_count
    errRate = error_count.astype(float) / y.shape[0]
    values = [sensitivity, specificity, PPV, NPV, errRate]
    names = ["sensitivity", "specificity", "PPV", "NPV", "errRate"]
    measures = zip(names, values)
    return measures

def get_performance_vals(y, pred_prob):
    """
        Use actual y and predicted probs to get performance measures
    """
    lowest_errRate = 1.
    #For measures that depend on threshold, choose threshold that minimizes errRate.  
    #Take lowest threshold when multiple thresholds give same minimum errRate.
        #Thus, we prioritize sensitivity over specificity, which is appropriate.
    for threshold in np.arange(0, 1.1, .1):
        #print "Threshold value: " , threshold
        performance_vals = measure_performance(y, pred_prob, threshold)
        #print performance_vals
        if performance_vals[-1][1] < lowest_errRate:
            #print "current error: " , performance_vals[-1][1]
            lowest_errRate = performance_vals[-1][1]
            best_performance = performance_vals
            best_threshold = threshold
    #meaures that are independent of threshold
    AUC = metrics.roc_auc_score(y, pred_prob)
    MSE = sum((y - pred_prob)**2)
    best_performance.append(("AUC", AUC))
    best_performance.append(("MSE", MSE))
    best_performance.append(("threshold", best_threshold))
    return best_performance

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
    #predictors = ["is.female", "Temperatura","Tos","Albumina"] #for testing
    #print "Predictors to use:\n" , predictors

    #create pandas dataframe with data that was cleaned in R
    df = pd.read_csv(inputsDir + "clin24_full_wImputedRF1.txt", sep='\t', 
        true_values=["True","Yes"], false_values=["False","No"], na_values=['NaN','NA']) 
    #df["is.DEN"] = df["is.DEN"].astype(int) #ML functions give predicted probs only for int outcomes
    df[["is.female"]] = (df[["Sexo"]]=="female") #need this to be True/False binary
    print "Number of rows in dataframe: " , df.shape[0], "\n"
    #print "Column names in dataframe: " , list(df.columns.values), "\n"  
    #import pdb 
    #pdb.set_trace()
    X = df[predictors].astype(float).values #include only vars we will use and convert to np array 
    y = df[outcome].astype(int).values #make outcome 0/1 and convert to np array
    #X, y=datasets.make_classification(n_samples=88, n_features=95) #generate toy data for testing

    #test functions with just 1 algorithm
    pred_prob = get_predictions_svm(X, y)
    #print zip(y, pred_prob)
    print "Performance results: " , get_performance_vals(y, pred_prob)    

    #run Super Learner
    RF = RandomForestClassifier()
    logitL1=linear_model.LogisticRegression(penalty='l1', solver='liblinear') 
        #todo: explore l1 path options for optimizing penalty coefficient (e.g., l1_min_c)
    logitL2=RidgeClassifier(normalize=True) #outputs class predictions but not probs
    logitL2.fit(X,y)
    pred_prob2 = logitL2.predict_proba(X)
    #print zip(pred_prob2, y)
    print "Performance results: " , get_performance_vals(y, pred_prob2)

    nn=neighbors.KNeighborsClassifier(n_neighbors=4) #Classes are ordered by lexicographic order.
    gradB_dev=GradientBoostingClassifier(loss='deviance') #for probabilistic outputs
    adaBoost=AdaBoostClassifier() 
    myLDA = LDA(solver='lsqr', shrinkage='auto') #optimal shrinkage is calculated analytically
    myQDA = QDA() 
    svm1=svm.SVC(kernel='rbf', probability=True)     
    #svm2=svm.SVC(kernel='poly', probability=True) #much too slow!
    svmL2=svm.LinearSVC(penalty='l2') #penalized linear support vector classification
    #todo: spectral clustering classifier?
    #todo: write a means classifier myself
    #lib=[RF, adaBoost]
    lib=[RF, logitL1, logitL2, nn, gradB_dev, adaBoost, myLDA, myQDA, svm1, svmL2]
    #libnames=["RF", "LogitL1", "LogitL2", "Nearest Neighbor","Gradient Boosting",
    #           "AdaBoost", "LDA", "QDA", "SVM-rbf", "SVM-L2"]
	
    #sl=SuperLearner(lib, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
    #sl.fit(X, y)
    #sl.summarize()
    #print "actual and predicted values"
    #print type(sl.y_pred_cv)
    #print sl.y_pred_cv.shape #dimension is n (#obs) by k (#algorithms)
    #print sl.y_pred_cv.shape
    start_time_cvSL = time.time()
    #cv_superlearner(sl, X, y, K=2, stratifyCV=True) #currently only returns risks_cv np array
    log_statement("cvSL execution time: {} minutes".format(
        (time.time() - start_time_cvSL)/60. ) ) 

    log_statement("Total execution time: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()