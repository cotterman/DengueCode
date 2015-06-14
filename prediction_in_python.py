###############################################################################
###### Employ machine-learning algorithms in Python to predict Dengue  ########
###############################################################################


import os  #has several functions for manipulating files and directories
import sys
import pickle
import warnings
import time
import pdb #debugger

import numpy as np
import pandas as pd
print "Version of pandas: " , pd.__version__ #should be v0.16.1
from pandas.core.categorical import Categorical

import matplotlib.pyplot as plt
import seaborn as sns #for prettier plots
from matplotlib import style
from matplotlib_style_utils import rstyle

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
import sklearn.cross_validation as cv
from sklearn.cross_validation import cross_val_score, cross_val_predict
from sklearn import grid_search, metrics, datasets
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn import svm, linear_model, neighbors
#from sklearn.linear_model import LogisticRegressionCV

sys.path.append('../SuPyLearner/supylearner')
import core #this is the main SuPyLearner code
from core import SuperLearner, cv_superlearner
from cross_val_utils import cross_val_predict_proba


np.random.seed(100)

#inputsDir = "/home/carolyn/dengue_data_and_results_local/intermediate_data/" #home PC
#inputsDir = "/home/nboley/Desktop/dengue_data_and_results_local/intermediate_data/" #home PC, diff login
inputsDir = "/srv/scratch/ccotter/intermediate_data/" #mitra

VERBOSE = True
def log_statement(statement):
    if VERBOSE: print statement


class RidgeClassifier(linear_model.RidgeClassifier):
    def predict_proba(self, X): 
        coefs = self.coef_
        z = np.dot(X, np.transpose(coefs))
        predict_proba1 = 1.0 / (1.0 + np.exp(-z)) #probability of getting 1
        predict_proba0 = 1.0 - predict_proba1     #probability of getting 0
        # return the predicted probabilities in same format as other learners do
        predict_proba = np.hstack([predict_proba0,predict_proba1])
        return predict_proba

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
    fitted_values = clf.predict_proba(X)
    #preds_cv = cross_val_predict(clf, X, y)
    print clf.classes_ #returns array containing order classes appear in fitted_values (todo)
    return fitted_values

def measure_performance(y, pred_prob, threshold, method_name=""):
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
    sensitivity = true_positives_count.astype(float) / positives_count
    specificity = true_negatives_count.astype(float) / negatives_count
    PPV = true_positives_count.astype(float) / predicted_positives_count 
    NPV = true_negatives_count.astype(float) / predicted_negatives_count
    errRate = error_count.astype(float) / y.shape[0]
    measures = pd.DataFrame({'sensitivity': [sensitivity], 
                             'specificity': [specificity],
                             'PPV': [PPV],
                             'NPV': [NPV],
                             'errRate': [errRate]}, 
                            index=[method_name])
    return measures

def get_cvAUC_with_CIs(y, pred_prob, cv_gen, confidence=0.95, method_name=""): 
    """
    Parameters
        y : array containing int 0s and 1s.  These are the true class labels.
        y_pred_prob: array containing predicted probabilities that y is 1. 
        cv_gen : cross-validation generator.  
             Should correspond to that which was used to generate y_pred_prob
        confidence: indicates desired coverage of confidence interval. 
                    Default = 0.95
    Returns
        list(cvAUC, se, ci, confidence)
    """    
    #Convert y array to array of ints if not already (so avoid rounding problems)
    if y.dtype=='int': y = y.astype(float)
    #create inverse probability weights to be used in ci.cvAUC calc
    n_obs = y.shape[0] #tot number of observations
    positives_count = (np.absolute(y-1)<1e-6).astype(int).sum()
    negatives_count = (np.absolute(y)<1e-6).astype(int).sum()
    w1 = 1.0 / (positives_count / float(n_obs))
    w0 = 1.0 / (negatives_count / float(n_obs)) 
    tracker = 0
    K = cv_gen.n_folds
    AUCs = np.zeros(K, dtype=float)
    ICvals = np.zeros(n_obs, dtype=float) #an array to hold influence curve values
    for train_index, test_index in cv_gen:
        y_test = y[test_index]
        pred_test = pred_prob[test_index]
        n_row = y_test.shape[0] #number of obs in current test set
        n_pos = (y_test==1).sum() #number of positive obs in test set       
        n_neg = n_row - n_pos #number of neg. obs in test set
        AUCs[tracker] = metrics.roc_auc_score(y_test, pred_test)
        for index in test_index:
            #when current y==1, count obs in test set for which y==0 and pred_y<pred_y_current
            if y[index]==1:
                fracNegYsWithSmallerPreds = np.logical_and(
                    y_test==0, pred_test < pred_prob[index]).astype(int).sum() / float(n_neg)
                ICvals[index] = w1 * (fracNegYsWithSmallerPreds - AUCs[tracker])
            #when current y==0, count obs in test set for which y==1 and pred_y>pred_y_current
            if y[index]==0:
                fracPosYsWithLargerPreds = np.logical_and(
                    y_test==1, pred_test > pred_prob[index]).astype(int).sum() / float(n_pos)
                ICvals[index] = w0 * (fracPosYsWithLargerPreds - AUCs[tracker])
        tracker += 1
    sighat2 = (ICvals**2).mean() / float(K)
    se = (sighat2 / n_obs)**.5
    cvAUC = AUCs.mean()
    z = stats.norm.ppf(.95)
    ci_low = max(cvAUC - (z*se), 0) #truncate CI at [0,1]
    ci_up = min(cvAUC + (z*se), 1) 
    cvAUCout = pd.DataFrame({'cvAUC': [cvAUC], 
                         'se': [se],
                         'ci_low': [ci_low],
                         'ci_up': [ci_up],
                         'confidence': [confidence]}, 
                        index=[method_name]) 
    return cvAUCout

def get_performance_vals(y, pred_prob, method_name, cv_gen, confidence=0.95):
    """
        Use actual y and predicted probs to get performance measures
    """
    lowest_errRate = 1.
    #For measures that depend on threshold, choose threshold that minimizes errRate.  
    #Take lowest threshold when multiple thresholds give same minimum errRate.
        #Thus, we prioritize sensitivity over specificity, which is appropriate.
    for threshold in np.arange(0, 1.1, .1):
        #print "Threshold value: " , threshold
        performance_vals = measure_performance(y, pred_prob, threshold, method_name)
        #print performance_vals['errRate']
        if float(performance_vals['errRate']) < lowest_errRate:
            #print "current error: " , performance_vals[-1][1]
            lowest_errRate = float(performance_vals['errRate'])
            best_performance = performance_vals
            best_threshold = threshold
    #meaures that are independent of threshold
    AUC = metrics.roc_auc_score(y, pred_prob)
    MSE = sum((y - pred_prob)**2)
    best_performance['AUC'] = AUC
    best_performance['MSE'] = MSE
    best_performance['threshold'] = best_threshold
    #cvAUC and confidence intervals
    cvAUCout = get_cvAUC_with_CIs(y, pred_prob, cv_gen, confidence, method_name)
    #combine results
    myperformance = pd.concat([best_performance, cvAUCout], axis=1)
    return myperformance

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
    #print "Actual outcomes: " , y[:10]
    #X, y=datasets.make_classification(n_samples=88, n_features=95) #generate toy data for testing

    if (True==False):
        #test functions with just 1 algorithm
        svm_preds = get_predictions_svm(X, y)[:, 0] #says class order is 0,1 but 0th ele is better 
        #print zip(y, svm_preds)
        cv_gen = cv.StratifiedKFold(y, n_folds=2, shuffle=True)
        svm_out = get_performance_vals(y, svm_preds, "mySVM", cv_gen, confidence=0.95)
        print "\nPerformance results using SVM:\n" , svm_out
        
        logitL2=RidgeClassifier(normalize=True) #outputs class predictions but not probs
        logitL2.fit(X,y)
        pred_prob = logitL2.predict_proba(X)[:, 1] #not cross-validated
        ridge_preds = cross_val_predict_proba(logitL2, X, y, cv_gen)[:, 1] 
        ridge_out = get_performance_vals(y, ridge_preds, "myRidge", cv_gen, confidence=0.95)
        print "\nPerformance results using Ridge: " , ridge_out

    start_time_cvSL = time.time()
    #run Super Learner
    cv_gen = cv.StratifiedKFold(y, n_folds=2, shuffle=True)
    RF = RandomForestClassifier()
    logitL1=linear_model.LogisticRegression(penalty='l1', solver='liblinear') 
        #todo: explore l1 path options for optimizing penalty coefficient (e.g., l1_min_c)
    logitL2=RidgeClassifier(normalize=True) #outputs class predictions but not probs
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
    lib=[RF, logitL1, logitL2, nn, gradB_dev, adaBoost, myLDA, myQDA, svm1]
    libnames=["RF", "Lasso", "Ridge", "NN","Gradient Boost",
               "AdaBoost", "LDA", "QDA", "SVM-rbf"]
    #test
    #for est in lib:
    #    est.fit(X, y)
    #    print est.classes_ #all look the same ( [0, 1] )
    sl=SuperLearner(lib, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
    SL_preds = cross_val_predict_proba(sl, X, y, cv_gen)
    resultsDF = get_performance_vals(y, SL_preds, "Super Learner", cv_gen, confidence=0.95)
    
    for counter, est in enumerate(lib):
        est.fit(X,y)
        preds = cross_val_predict_proba(est, X, y, cv_gen)[:, 0]
        out = get_performance_vals(y, preds, libnames[counter], cv_gen, confidence=0.95)
        #merge results together
        resultsDF = pd.concat([resultsDF, out])
    print "\nCombined results: \n" , resultsDF
    #sort results by decreasing cvAUC
    resultsDF.sort(columns=['cvAUC','errRate'], axis=0, ascending=False, inplace=True)

    #print "SL cv predictions: " , SL_preds[:10]
    #print "\nPerformance results, SL: " , get_performance_vals(y, SL_preds)
    #print sl.y_pred_cv.shape # numpy array of dimension n (#obs) by k (#algorithms)
    #print "\nPerformance results, RF: " , get_performance_vals(y, sl.y_pred_cv[:,0])

    #make plots
    labels = resultsDF.index.values #array with method labels
    label_pos = np.arange(len(labels))
    measurements = np.array(resultsDF['cvAUC'])
    #plotting functions want se not confidence interval bounds
    z = stats.norm.ppf(.95)
    SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
           ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]

    #use matplotlib defaults
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.barh(label_pos, measurements, xerr=SEs, align='center', alpha=0.4)
    plt.yticks(label_pos, labels)
    plt.xlabel('cvAUC')
    plt.title('Distinguishing OFI from DENV')
    style.use('ggplot') #alternative 1 -- use defaults that mimic R's ggplot
    rstyle(ax) #alternative 2 -- gives same result as style.use('ggplot')
    plt.show()

    #alternative 2 -- use seaborn
    #problem: it seems to want to calculate error bars and not allow me to give them
    #sns.set(style="white", context="talk")
    #sns.set(style='ticks', palette='Set2') #prettyplot recommendation
    #f, ax1 = plt.subplots(1, 1, figsize=(8, 6))
    #sns.barplot(labels, measurements, ci=SEs, ax=ax1)
    #ax1.set_ylabel("cvAUC")
    #sns.despine()
    #plt.setp(f.axes, yticks=[])
    #plt.show()

    #cv_superlearner(sl, X, y, K=2, stratifyCV=True) #currently only returns risks_cv np array
    log_statement("\ncvSL execution time: {} minutes".format(
        (time.time() - start_time_cvSL)/60. ) ) 

    log_statement("Total execution time: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()
