###############################################################################
###### Employ machine-learning algorithms in Python to predict Dengue  ########
###############################################################################


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

inputsDir = "/home/ccotter/dengue_data_and_results_local/intermediate_data/" #home PC
outDir = "/home/ccotter/dengue_data_and_results_local/python_out/" #home PC
#inputsDir = "/home/nboley/Desktop/dengue_data_and_results_local/intermediate_data/" #home PC diff login
#inputsDir = "/srv/scratch/ccotter/intermediate_data/" #mitra
#outDir = "/srv/scratch/ccotter/python_out/" #mitra

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

def get_predictor_desc(file_name, outcome):
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
    if outcome == "is.DHF_DSS" and predictors.count("IR")>0 :
        predictors.remove("IR") #todo: turn to binary int variable
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

def get_performance_vals(y, pred_prob, method_name, cv_gen, confidence, startingDF):
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
    #merge results together with the DF that was passed
    resultsDF = pd.concat([startingDF, myperformance])
    return resultsDF

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

def get_data(inputsDir, filename, NoInitialDHF, patient_sample):
    df = pd.read_csv(inputsDir + filename, sep='\t', 
        true_values=["True","Yes"], false_values=["False","No"], na_values=['NaN','NA']) 
    #df["is.DEN"] = df["is.DEN"].astype(int) #ML funcs give predicted probs only for int outcomes
    df[["is.female"]] = (df[["Sexo"]]=="female").astype(int)
    df[["is.cohort"]] = (df[["Study"]]=="Cohort").astype(int) 
    if (NoInitialDHF==True):
        df = df[df.WHO12hr4cat!="DSS"] 
        df = df[df.WHO12hr4cat!="DHF"] 
        df = df[df.DENV=="Positivo"] #limit to only DENV positive patients
    if (patient_sample == "hospital_only"):
        df = df[df.Study=="Hospital"] #limit to only hospital patients
    if (patient_sample == "cohort_only"):
        df = df[df.Study=="Cohort"] #limit to only cohort patients
    print "Number of rows in dataframe: " , df.shape[0], "\n"
    #print "Column names in dataframe: " , list(df.columns.values), "\n" 
    return df

def add_imput_dummies(include_imp_dums, df, predictors):
    """
        Find imputation dummies in data that correspond to predictor list
        Add these to the list of predictors if not redundant with is.cohort indicator
            If redundant with is.cohort, ensure that is.cohort is in data  
    """
    if (include_imp_dums==True):
        allvars = df.columns.values #np.array
        imp_found = []
        for var in allvars:
            y = re.findall(".+?_imputed", var)
            assert len(y) in (0,1)
            if len(y)==1: imp_found.append(y[0])
        imp_imagine = [var+'_imputed' for var in predictors] 
        imp_matched = [var for var in imp_found if var in imp_imagine]   
        #if there are imput vars that are redundant with is.cohort, drop them
            # and add is.cohort if not already included.  this will improve interpretation of results
        cohort_imitator = 0
        for var in imp_matched:
            if sum((df[[var]].values==df[["is.cohort"]].values)==False)[0] == 0:
                imp_matched.remove(var) #remove if perfect match
                cohort_imitator = 1 #keep track of fact that we found a match
        if cohort_imitator==1:
            if 'is.cohort' not in imp_matched: imp_matched.append('is.cohort')
        predictors = predictors + imp_matched    
    return predictors

def build_library(p, nobs, screen=None):
    """
        Develop list of prediction functions.
        p is the number of predictors.  
        nobs in number of observations in data
            The parameters fed to some algorithms (e.g., RandomForest) are functions of p and nobs.
    """
    print "Number of vars: " , p

    #set parameters for GridSearchCV
    n_jobs = 5
    cv_grid = 5

    mean = dummy.DummyClassifier(strategy='most_frequent') #predict most frequent class

    #very simple CART (for comparison sake)
    CART = tree.DecisionTreeClassifier(max_depth=10) 

    #todo: consider adding a shrinkage threshold (chosen via gridsearch)
    centroids = neighbors.NearestCentroid(metric='euclidean', shrink_threshold=None)

    LDA_shrink = LDA(solver='lsqr', shrinkage='auto') #optimal shrinkage is calculated analytically

    myQDA = QDA()

    GBparams = {'max_depth':[2,3,4], 'min_samples_split':[2,4], 'min_samples_leaf':[1,3],
        'max_features':[None], 'max_leaf_nodes':[None], 'loss':['deviance']}
    GBtune = grid_search.GridSearchCV(GradientBoostingClassifier(), GBparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    adaBoost=AdaBoostClassifier() 

    #n_estimators is number of trees in forest (R default is 1000)
    #max_features is the number of features to consider when looking for best split
    if p < 9:
        RF_max_features = [ int(math.ceil(math.sqrt(p))) ]
    else:
        RF_max_features = [ int(math.ceil(math.sqrt(p))), int(math.ceil(p/3)) ]
    RFparams = {'n_estimators':[1000],  'max_features': RF_max_features}
 
    RFtune = grid_search.GridSearchCV(RandomForestClassifier(), RFparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    NNparams = {'n_neighbors':[3,5,7]}
    NNtune = grid_search.GridSearchCV(neighbors.KNeighborsClassifier(), NNparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    #C is the inverse of regularization strength.  Default is 1.
        #alpha (parameter in SGD and Ridge) is defined as (2*C)^-1
    L1params = {'penalty':['l1'], 'C':[.5, 1, 1.5], 'solver':['liblinear']}
    L1tune = grid_search.GridSearchCV(linear_model.LogisticRegression(), L1params,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid) 

    #larger alpha --> more regularization (larger penalty on complexity term)
    #In practice, performance did not differ with standardized data
    if nobs+20 > p:
        alpha_start = .0001
    else:
        alpha_start = .001 #need greater dimension reduction when p>>n 
    L2params = {'alpha':list(np.linspace(start=alpha_start, stop=2, num=100))}
    L2tune = grid_search.GridSearchCV(RidgeClassifier(normalize=True), L2params,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid) 

    #l1_ratio of 1 corresponds to L1 penalty (lasso); l1_ratio of 0 gives L2 penalty (ridge)
        #note: R's the elastic net penalty is [(1-param)/2]L2 + paramL1 
        #Rather than simply (1-param)L2 + paramL1 as it is here (Tibshironi does as python does) 
    #alpha is the regularization term coefficient (called lambda in R's glmnet)
    #Does better with standardized data
    if nobs+20 > p:
        alpha_start = .0001
    else:
        alpha_start = .001 #need greater dimension reduction when p>>n 
    ENparams = {'l1_ratio':[0, .15, .3, .5, .7, .85, 1], 'loss':['log'], 'penalty':['elasticnet'],
                'alpha':list(np.linspace(start=alpha_start, stop=2, num=100)), 
                'warm_start':[True], 'fit_intercept':[True]} 
    ENtune = grid_search.GridSearchCV(linear_model.SGDClassifier(), ENparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    #todo: try ensemble gradient descent

    #uses squared hinge loss by default
    #C is the inverse regularization path
    SVMparams = {'penalty':['l2'], 'C':[.5, 1, 1.5], 'fit_intercept':[True]}
    svmL2tune = grid_search.GridSearchCV(svm.LinearSVC(), SVMparams,
         score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid) 

    #hinge loss and rbf(radial basis function)
    RBFparams = {'kernel':['rbf'], 'C':[.5, 1, 1.5], 'probability':[True]}
    svmRBFtune = grid_search.GridSearchCV(svm.SVC(), RBFparams,
            score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid) 
     
    #libs=[CART, adaBoost, NNtune] #for testing
    #libnames = []
    libs=[mean, CART, centroids, LDA_shrink, myQDA, GBtune, adaBoost, RFtune, 
        NNtune,L1tune, L2tune, ENtune, svmL2tune, svmRBFtune]
    libnames=["mean", "CART", "centroids", "LDA.shrink", "QDA", "G.Boost", "adaBoost", "RF", 
        "NN", "logitL1", "logitL2", "E.Net", "svmL2", "svmRBF"]
    if libnames == []:
        libnames=[est.__class__.__name__ for est in libs]

    #add feature screen if specified
    if screen != None:
        for counter, lib in enumerate(libs):
            lib_screened = Pipeline([ ('feature_selection', screen), (libnames[counter], lib) ])
            libs[counter] = lib_screened

    return (libs, libnames)

def get_screen(screenType, screenNum, predictors):
    """
        Configure the screen (variable selection method) 
        for specified screenType and screenNum
    """
    if (screenType==None):
        screenNum = len(predictors)
        screen = None
    if (screenType=="univariate"):
        myscreen = SelectKBest(feature_selection.f_classif, k=screenNum) #works
    elif (screenType=="L1"):
        L1_filter = svm.LinearSVC(C=.01, penalty='l1', dual=False) #works #small C -> few selected
        #todo: output avg number of features selected using SVC
    elif (screenType=="RFE"):    
        #RFE first trains estimator using all features and assigns weights to each of them 
            #(e.g., weights may be the coefficients of a linear model)
            #features with smallest weights are pruned,
            #process repeated until desired number of of covars is obtained 
        #screen = feature_selection.RFE(screenNum) #recursive feature elimination -- todo
        print "RFE: todo"
    elif (screenType=="tree"):
        #screen = ExtraTreesClassifier().fit(X,y).transform(X).feature_importances_ -- todo
        print "tree screen: todo"
    elif (screenType=="L1random"):
        #perturb design matrix or sub-sample and count number of times given regressor is selected
            #mitigates the problem of selecting 1 feature from a group of very correlated ones
        #screen = RandomizedLogisticRegression -- todo
        print "randomL1 screen: todo"
    return [screen, screenNum]

def results_for_library(X, y, cv_gen, libs, libnames, resultsDF): 
    """
        Fit each algorithm in library to data and add performance measures
        to pandas dataframe (resultsDF)
    """    
    for counter, est in enumerate(libs): 
        print "Now fitting " , libnames[counter]       
        est.fit(X,y)
        if hasattr(est, "best_estimator_"):
            print "Optimal parameters: " , est.best_estimator_
        if hasattr(est, "predict_proba"):
            #exceptions for the algorithms that give results in diff order
            if est.__class__.__name__ == "SVC":
                preds = cross_val_predict_proba(est, X, y, cv_gen)[:, 0]
            else:
                preds = cross_val_predict_proba(est, X, y, cv_gen)[:, 1]
        else:
            preds = cross_val_predict(est, X, y, cv_gen)
        resultsDF = get_performance_vals(y, preds, libnames[counter], 
                                    cv_gen, 0.95, resultsDF)
    return resultsDF

def add_info_and_print(resultsDF, include_imp_dums, screenType, pred_count, 
                        patient_sample, nobs, FileNamePrefix, predictor_desc, outDir, printme):
    resultsDF.insert(loc=len(resultsDF.columns), column='include_imp_dums', value=include_imp_dums)
    resultsDF.insert(loc=len(resultsDF.columns), column='screen_method', value=screenType)
    resultsDF.insert(loc=len(resultsDF.columns), column='predictor_count', value=pred_count)
    resultsDF.insert(loc=len(resultsDF.columns), column='patient_sample', value=patient_sample)
    resultsDF.insert(loc=len(resultsDF.columns), column='n', value=nobs)
    #sort results by decreasing cvAUC
    resultsDF.sort(columns=['cvAUC','errRate'], axis=0, ascending=False, inplace=True)
    print "\nCombined results: \n" , resultsDF
    #output results to excel file
    if printme==True:
        outFileName = FileNamePrefix + '_' + predictor_desc + '.txt'
        resultsDF.to_csv(outDir+outFileName, sep=",")    
    #alternative: could use the predicted values coming from running SL
    #print sl.y_pred_cv.shape # numpy array of dimension n (#obs) by k (#algorithms)
    #print "\nPerformance results, RF: " , get_performance_vals(y, sl.y_pred_cv[:,0])
    return resultsDF

def plot_results(resultsDF, plot_title, figName, outDir):
    """
        Create plots of cvAUC along with error bars
    """
    labels = resultsDF.index.values #array with method labels
    label_pos = np.arange(len(labels))
    measurements = np.array(resultsDF['cvAUC'])
    #plotting functions want se not confidence interval bounds
    z = stats.norm.ppf(.95)
    SEs = [( np.array(resultsDF['cvAUC']) - np.array(resultsDF['ci_low']) )/z, 
           ( np.array(resultsDF['ci_up']) - np.array(resultsDF['cvAUC']) )/z ]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.barh(label_pos, measurements, xerr=SEs, align='center', alpha=0.4)
    plt.yticks(label_pos, labels)
    plt.xlabel('cvAUC')
    plt.title(plot_title)
    #style.use('ggplot') #alternative 1 -- use defaults that mimic R's ggplot
    #rstyle(ax) #alternative 2 -- gives same result as style.use('ggplot')
    plt.savefig(outDir + figName)
    plt.show() 

def main():
    start_time_overall = time.time()

    ## Choose outcome variable ##
    outcome = "is.DEN"  
    #outcome = "is.DHF_DSS"
    if outcome=="is.DEN":
        comparison_groups = "OFI vs. DENV using " #will appear in graph title
        FileNamePrefix = "OFI.v.DEN" #use all samples; exclude IR and PCR predictors
        NoInitialDHF = False #whether to exclude samples with initial DHF/DSS diagnosis
    elif outcome=="is.DHF_DSS":
        comparison_groups = "DF vs. DHF/DSS using " #will appear in graph title
        FileNamePrefix = "Predict_DF.v.DHFandDSS" #exclude OFI+early DHF/DSS; use all predictors
        NoInitialDHF = True #whether to exclude samples with initial DHF/DSS diagnosis 
    if (NoInitialDHF==True):
        restrictions = ", DENV patients with non-severe initial Dx"
    else:
        restrictions = ""

    ## Choose patient sample ##
    #patient_sample = "all"
    #patient_sample = "cohort_only"
    patient_sample = "hospital_only"
    
    ## Choose list of variables to use in prediction ##
    predictor_desc = "covarlist_all"
    #predictor_desc = "covarlist_noUltraX"
    #predictor_desc = "covarlist_CohortRestrict"
    #predictor_desc = "covarlist_genOnly"
    #predictors = ["is.female", "Temperatura","Tos","Albumina"] #for testing
    predictors = get_predictor_desc(predictor_desc+".txt", outcome)
    print "Original predictor list:\n" , predictors

    ## Create pandas dataframe with data that was cleaned in R ##
    df = get_data(inputsDir, "clin12_full_wImputedRF1.txt", NoInitialDHF, patient_sample)

    ## Choose whether to include indicator of clinic type (hospital v cohort) ##
    include_study_dum = True
    if patient_sample=="cohort_only" or patient_sample=="hospital_only":
        include_study_dum = False
    if (include_study_dum==True):
        predictors.append("is.cohort")

    ## Choose whether to include imputation dummies ##
    include_imp_dums = False
    #modify predictor list so as to include imputation dummies if desired
    predictors = add_imput_dummies(include_imp_dums, df, predictors)
    print "Predictors to include, pre-screening:\n" , predictors
 
    ## Choose variable screening method (if any) ##
    screenType = None #current options: None, "univariate", "L1"
    screenNum = 5 #applies when screenType != None; else will be ignored
    screen, pred_count = get_screen(screenType, screenNum, predictors)

    ## Build library of classifiers ##
    myLibrary = build_library( p=len(predictors), nobs=df.shape[0], screen=screenType)
    libs = myLibrary[0]
    libnames = myLibrary[1]

    ## Keep only columns in predictors list, create arrays of scaled variables ##
    X_prelim = df[predictors].astype(float).values #include only vars we will use and convert to np array
    #scale data since ridge/elastic net are not equivariant under scaling
        #other algorithms should be unaffected by scaling
        #In practice, standardization appears to work better than MinMaxScaler
        #X = preprocessing.MinMaxScaler(feature_range=(0,1)).fit_transform(X) #scale features to a range
    X = preprocessing.scale(X_prelim) #standardize data (subtract mean, then divide by std) 
    y = df[outcome].astype(int).values #make outcome 0/1 and convert to np array
    #print "Actual outcomes: " , y[:10]
    #X, y=datasets.make_classification(n_samples=88, n_features=95) #generate toy data for testing

    ## Get Predictions ##
    cv_gen = cv.StratifiedKFold(y, n_folds=2, shuffle=True, random_state=10)
    resultsDF = pd.DataFrame() #empty df to hold results
    # Super Learner
    start_time_cvSL = time.time()
    #sl = SuperLearner(libs, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
    #SL_preds = cross_val_predict_proba(sl, X, y, cv_gen)
    #resultsDF = get_performance_vals(y, SL_preds, "Super Learner", cv_gen, 0.95, resultsDF)
    end_time_cvSL = time.time()
    # Results for each algorith in library
    resultsDF = results_for_library(X, y, cv_gen, libs, libnames, resultsDF) 

    ## Add columns with additional methods info, print results to text file ##
    resultsDF = add_info_and_print(resultsDF, include_imp_dums, screenType, pred_count, 
                        patient_sample, df.shape[0], 
                        FileNamePrefix, predictor_desc, outDir, printme=True)

    ## Make plots of results ##
    plot_title = comparison_groups+predictor_desc+restrictions
    figName = FileNamePrefix + '_' + predictor_desc + '.png'
    plot_results(resultsDF, plot_title, figName, outDir)
    
    ## Ouput execution time info ##
    log_statement("\ncvSL execution time: {} minutes".format(
        (end_time_cvSL - start_time_cvSL)/60. ) ) 

    log_statement("Total execution time: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()
