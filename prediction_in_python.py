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

import matplotlib as mpl
mpl.use('Agg') #avoids running an X server (to avoid errors with remote runs)
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt; plt.rcdefaults()
#from matplotlib import style
#from matplotlib_style_utils import rstyle
import matplotlib.lines as mlines
import seaborn as sns #this will change many graphing parameters for style


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

def get_predictor_desc(file_name, outcome, NoOFI):
    """
        Imports list of variables contained in file_name
        and replaces categorical variable names with binary equivalents
    """
    if file_name=="covarlist_custom.txt":
        f = open(inputsDir + 'covarlist_genOnly.txt', 'r')
        #add top 4 lab vars
        predictors = f.read().splitlines() + ['Plaquetas','Leucocitos','TGP_','TGO_'] 
    else:   
        f = open(inputsDir + file_name, 'r')
        predictors = f.read().splitlines()

    #drop variables that are not appropriate for prediction type
    if outcome=="is.DEN" or NoOFI==False:
        if predictors.count("PCR")>0 :
            predictors.remove("PCR")
    #drop variables that are not appropriate for prediction type
    #if outcome=="is.DHF_DSS":
        
    #use binary equivalents of the categorical predictors
        #Torniquete -- less than 20 vs. 20+ = is.torniquete20plus
        #Pulso  -- strong or moderate vs. rapid or not palpable = is.pulse_danger
        #PCR -- 1, 2, or 3 serotype (DENV positive only) = is.serotype1,  etc.
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

    #should never include these variables in prediction analysis 
        #since they are obtained from convelescent samples
    if predictors.count("IR")>0 :
        predictors.remove("IR") 
    if predictors.count("HematocritoElev")>0:
        predictors.remove("HematocritoElev") 
    if predictors.count("Hemoconcentracion")>0:
        predictors.remove("Hemoconcentracion")        
        
    return(predictors)


def results_for_testset(X, y, Xnew, ynew, libs, libnames, sl_folds=5): 
    """
        Fit each algorithm in library to (X, y) and predict for Xnew
        add preds and performance measures to pandas dataframes 
        (predDFnew and resultsDFnew)
    """    
    predDFnew = pd.DataFrame() #this will hold predicted probs for all algorithms
    resultsDFnew = pd.DataFrame() #empty df to hold performance measures
    # add super learner to library list
    sl = SuperLearner(libs, loss="nloglik", K=sl_folds, stratifyCV=True, save_pred_cv=True)
    libs.append(sl)
    libnames.append('Super Learner')
    # fit (to hospit data) and predict (on cohort data)
    for counter, est in enumerate(libs): 
        print "Now fitting " , libnames[counter]       
        est.fit(X,y)
        if hasattr(est, "best_estimator_"):
            print "Optimal parameters: " , est.best_estimator_
        if libnames[counter]=='Super Learner':
            preds = est.predict_proba(Xnew)
        elif hasattr(est, "predict_proba"):
            #exceptions for the algorithms that give results in diff order
            if est.__class__.__name__ == "SVC":
                preds = est.predict_proba(Xnew)[:, 0]
            else:
                preds = est.predict_proba(Xnew)[:, 1]
        else:
            preds = est.predict(Xnew)
        #save the predicted values for each algorithm in one dataframe
        predDFnew.insert(loc=len(predDFnew.columns), column=libnames[counter], value=preds)
        #get performance measures (add row to resultsDF)
        resultsDFnew = get_performance_vals(ynew, preds, 
                        libnames[counter], None, None, resultsDFnew)
    return predDFnew, resultsDFnew

def measure_performance(y, pred_prob, threshold, method_name=""):
    """
        Obtain measures that depend on specified classification threshold
    """
    pred_y = (pred_prob>threshold).astype(int) #turn predicted probs into 0/1 value
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
            #when current y==1, count obs in test set 
                #for which y==0 and pred_y<pred_y_current
            if y[index]==1:
                fracNegYsWithSmallerPreds = np.logical_and(
                    y_test==0, 
                    pred_test < pred_prob[index]).astype(int).sum() / float(n_neg)
                ICvals[index] = w1 * (fracNegYsWithSmallerPreds - AUCs[tracker])
            #when current y==0, count obs in test set 
                #for which y==1 and pred_y>pred_y_current
            if y[index]==0:
                fracPosYsWithLargerPreds = np.logical_and(
                    y_test==1, 
                    pred_test > pred_prob[index]).astype(int).sum() / float(n_pos)
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

def get_performance_vals(y, pred_prob, method_name, cv_gen, 
                        confidence, startingDF):
    """
        Use actual y and predicted probs to get performance measures
    """
    lowest_errRate = 1.
    #For measures that depend on threshold, choose threshold that minimizes errRate.  
    #Take lowest threshold when multiple thresholds give same minimum errRate.
        #Thus, we prioritize sensitivity over specificity, which is appropriate.
    for threshold in np.arange(0, 1.01, .01):
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
    #no need for cvAUC when we are finding predictions for test set
    if cv_gen != None:
        #cvAUC and confidence intervals
        cvAUCout = get_cvAUC_with_CIs(y, pred_prob, cv_gen, confidence, method_name)
        #combine results
        myperformance = pd.concat([best_performance, cvAUCout], axis=1)
        #merge results together with the DF that was passed
        resultsDF = pd.concat([startingDF, myperformance])
    else:
        resultsDF = pd.concat([startingDF, best_performance])
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

def add_imput_dummies(include_imp_dums, imp_dums_only, df, predictors_prelim):
    """
        Find imputation dummies in data that correspond to predictor_prelim list
        Add these to the list of predictors_prelim if not redundant with is.cohort indicator
            If redundant with is.cohort, ensure that is.cohort is in data  
        If imp_dums_only is true, then return predictor list containing imputation dummies only
    """
    if (include_imp_dums==True):
        allvars = df.columns.values #np.array
        imp_found = []
        for var in allvars:
            y = re.findall(".+?_imputed", var)
            assert len(y) in (0,1)
            if len(y)==1: imp_found.append(y[0])
        imp_imagine = [var+'_imputed' for var in predictors_prelim] 
        imp_matched = [var for var in imp_found if var in imp_imagine]   
        #if there are imput vars that are redundant with is.cohort, drop them
            # and add is.cohort if not already included.  
            # this may improve interpretation of results for things like RF
        cohort_imitator = 0
        for var in imp_matched:
            if sum((df[[var]].values==df[["is.cohort"]].values)==False)[0] == 0:
                imp_matched.remove(var) #remove if perfect match
                cohort_imitator = 1 #keep track of fact that we found a match
        if cohort_imitator==1:
            if 'is.cohort' not in imp_matched: imp_matched.append('is.cohort')
        if imp_dums_only==True:
            predictors = imp_matched
        else:
            predictors = predictors_prelim + imp_matched 
    else:
        predictors = predictors_prelim

    return predictors

def get_data(inputsDir, filename, NoInitialDHF, patient_sample, 
            NoOFI, outcome, predictors_prelim, include_study_dum, 
            include_imp_dums, imp_dums_only, standardize=True):
    """
        Create pandas dataframe from text file created in R
        Standardize predictor variables (if standardize=True)
        Eliminate columns according to:
            predictors_prelim (List of covariates, excluding study and imputation indictors)
            include_study_dum (True to include is.cohort variable)
            include_imp_dums (True to include imputation dummy variables)
        Eliminate observations according to:
            NoInitialDHF: True means to exclude patients with initial severe dengue Dx
            patient_sample: "cohort_only","hospital_only" or "both"
            NoOFI: True means to eliminate non-dengue patients
    """
    df_prelim = pd.read_csv(inputsDir + filename, sep='\t', 
        true_values=["True","Yes","TRUE"], false_values=["False","No","FALSE"], 
        na_values=['NaN','NA']) 
    #ML funcs give predicted probs only for int outcomes
        #currently turning to int as part of read_csv
    #df["is.DEN"] = df["is.DEN"].astype(int) 
    df_prelim[["is.female"]] = (df_prelim[["Sexo"]]=="female").astype(int)
    df_prelim[["is.cohort"]] = (df_prelim[["Study"]]=="Cohort").astype(int) 

    #standardize data before removing rows 
        #want standardization to be same for cohort and hospit data 
        #scale data since ridge/elastic net are not equivariant under scaling
    X_prelim = df_prelim[predictors_prelim]
    if (standardize==True):
        X = X_prelim.apply(lambda x: (x - x.mean()) / x.std() )
    else:
        X = X_prelim

    #modify predictor list so as to include cohort and imputation dummies if desired
        #note we are adding them after standardizing to prevent them from being standardized
    if include_study_dum==True:
        predictors = predictors_prelim + ["is.cohort"]
    else:
        predictors = predictors_prelim
    all_predictors = add_imput_dummies(include_imp_dums, imp_dums_only, df_prelim, predictors)
    #get list of variables that were not scaled (could be an empty list)
    addons =  [v for v in all_predictors if v not in predictors_prelim]  

    #put back the variables that were not scaled
    nonXdf = df_prelim[['code','Study','DENV','WHO_initial_given',outcome] + addons]
    if imp_dums_only == True: #only use imputation indicators
        df = nonXdf 
    else: #include clinical info plus possible addons
        df = pd.concat([nonXdf, X], axis=1)
    
    #remove rows according to parameter values
    if (NoInitialDHF==True):
        df = df[df.WHO_initial_given!="DSS"] 
        df = df[df.WHO_initial_given!="DHF"] 
    if (NoOFI==True):
        df = df[df.DENV=="Positivo"] #limit to only DENV positive patients
    if (patient_sample == "hospital_only"):
        df = df[df.Study=="Hospital"] #limit to only hospital patients
    if (patient_sample == "cohort_only"):
        df = df[df.Study=="Cohort"] #limit to only cohort patients
    print "Number of rows in dataframe: " , df.shape[0], "\n"
    #print "Column names in dataframe: " , list(df.columns.values), "\n" 

    return df, all_predictors

def build_library(p, nobs, screen=None, testing=False, univariate=False):
    """
        Develop list of prediction functions.
        p is the number of predictors.  
        nobs is number of observations in data
            Parameters fed to some algorithms (e.g., RF) are functions of p and nobs.
    """
    print "Number of vars: " , p

    #set parameters for GridSearchCV
    n_jobs = 5
    cv_grid = 5

    mean = dummy.DummyClassifier(strategy='most_frequent') #predict most frequent class

    #very simple CART (for comparison sake)
    CARTparams = {'max_depth': list(np.arange(21)[1:]) } #try max_depths of 1 thru 20
    CARTune = grid_search.GridSearchCV(tree.DecisionTreeClassifier(), CARTparams)

    #todo: consider adding a shrinkage threshold (chosen via gridsearch)
    centroids = neighbors.NearestCentroid(metric='euclidean', shrink_threshold=None)

    #optimal shrinkage is calculated analytically
    LDA_shrink = LDA(solver='lsqr', shrinkage='auto') 

    myQDA = QDA() #produces terrible results - todo: figure out why and add back in

    GBparams = {'max_depth':[2,3,4], 'min_samples_split':[2,4], 
        'min_samples_leaf':[1,3], 'n_estimators':[100],
        'max_features':[None], 'max_leaf_nodes':[None], 'loss':['deviance']}
    GBtune = grid_search.GridSearchCV(GradientBoostingClassifier(), GBparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    adaBoost=AdaBoostClassifier(n_estimators=500) 

    #n_estimators is number of trees in forest (R default is 1000)
    #max_features is the number of features to consider when looking for best split
    if p < 9:
        RF_max_features = [ int(math.ceil(math.sqrt(p))) ]
    else:
        RF_max_features = [ int(math.ceil(math.sqrt(p))), int(math.ceil(p/3)) ]
    RFparams = {'n_estimators':[500],  'max_features': RF_max_features}
    RFtune = grid_search.GridSearchCV(RandomForestClassifier(), RFparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    NNparams = {'n_neighbors':[3,5,7,9,11,13], 'weights':['uniform','distance']}
    NNtune = grid_search.GridSearchCV(neighbors.KNeighborsClassifier(), NNparams,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

    #C is the inverse of regularization strength.  Default is 1.
        #alpha (parameter in SGD and Ridge) is defined as (2*C)^-1
    L1params = {'penalty':['l1'], 'C':[.5, 1, 1.3, 1.5, 1.7], 'solver':['liblinear']}
    L1tune = grid_search.GridSearchCV(linear_model.LogisticRegression(), L1params,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid) 

    #larger alpha --> more regularization (larger penalty on complexity term)
    #In practice, performance did not differ with standardized data
    if nobs+20 > p:
        alpha_start = .0001
    else:
        alpha_start = .0001 #need greater dimension reduction when p>>n 
    L2params = {'alpha':list(np.linspace(start=alpha_start, stop=1, num=100))}
    L2tune = grid_search.GridSearchCV(RidgeClassifier(normalize=True), L2params,
        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid) 

    #l1_ratio of 1 corresponds to L1 penalty (lasso)
    #l1_ratio of 0 gives L2 penalty (ridge)
        #note: R's elastic net penalty is [(1-param)/2]L2 + paramL1 
        #Python and Tibshironi use (1-param)L2 + paramL1
    #alpha is the regularization term coefficient (called lambda in R's glmnet)
    #Does better with standardized data
    if nobs+20 > p:
        alpha_start = .0001
    else:
        alpha_start = .001 #need greater dimension reduction when p>>n 
    ENparams = {'l1_ratio':[0, .15, .3, .5, .7, .85, 1], 'loss':['log'], 
                'penalty':['elasticnet'],
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

    #algorithms that only apply to small p problems (e.g., for univariate VIM)
    logistic = linear_model.LogisticRegression(penalty=None) #gave errors
    LDA_noshrink = LDA(shrinkage=None)   
     
    if (testing==True):
        libs=[tree.DecisionTreeClassifier()] #for testing
        libnames = ["tree"]
    elif (univariate==True):
        #CART = tree.DecisionTreeClassifier()
        #libs=[CART, CART]
        #libnames=["C1","C2"]
        libs=[CARTune, mean, LDA_noshrink, myQDA, NNtune]
        libnames=["CART","Mean","LDA","QDA","N.Neighbor"] 
    else:
        libs=[mean, CARTune, centroids, LDA_shrink, GBtune, adaBoost, RFtune, 
            NNtune, L1tune, L2tune, svmL2tune] #svmRBFtune, ENtune -- wacky results
        libnames=["Mean", "CART", "Centroids", "LDA+shrinkage", "Gradient Boost", 
            "AdaBoost", "Random Forests", 
            "N.Neighbor", "Logit-L1", "Logit-L2", "SVM-L2"] #"SVM-RBF", "Elastic Net"
    if libnames == []:
        libnames=[est.__class__.__name__ for est in libs]

    #add feature screen if specified
    if screen != None:
        for counter, lib in enumerate(libs):
            lib_screened = Pipeline([ ('feature_selection', screen), 
                                        (libnames[counter], lib) ])
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
        #small C -> few selected
        L1_filter = svm.LinearSVC(C=.01, penalty='l1', dual=False) #works 
        #todo: output avg number of features selected using SVC
    elif (screenType=="RFE"):    
        #RFE first trains estimator using all features and assigns weights to each 
            #(e.g., weights may be the coefficients of a linear model)
            #features with smallest weights are pruned,
            #process repeated until desired number of of covars is obtained 
        #screen = feature_selection.RFE(screenNum) # -- todo
        print "RFE: todo"
    elif (screenType=="tree"):
        #todo:
            #screen = ExtraTreesClassifier().fit(X,y).transform(X).feature_importances_
        print "tree screen: todo"
    elif (screenType=="L1random"):
        #perturb design matrix or sub-sample and 
            #count number of times regressor is selected
            #mitigates problem of selecting 1 feature from group of very correlated ones
        #screen = RandomizedLogisticRegression -- todo
        print "randomL1 screen: todo"
    return [screen, screenNum]

def results_for_library(X, y, cv_gen, libs, libnames, predDF, resultsDF): 
    """
        Fit each algorithm in library to data and add predicted probabilities
        and performance measures to pandas dataframes (predDF and resultsDF)
    """    
    for counter, est in enumerate(libs): 
        print "Now running CV for " , libnames[counter]       
        if hasattr(est, "predict_proba"):
            #exceptions for the algorithms that give results in diff order
            if est.__class__.__name__ == "SVC":
                preds = cross_val_predict_proba(est, X, y, cv_gen)[:, 0]
            else:
                preds = cross_val_predict_proba(est, X, y, cv_gen)[:, 1]
        else:
            preds = cross_val_predict(est, X, y, cv_gen)
  
        #optional: take a look at parameters chosen via grid search 
        est.fit(X,y)
        if hasattr(est, "best_estimator_"):
            print "Optimal estimator: " , est.best_estimator_
            print "Parameters chosen (via CV): " , est.best_params_ 

        #save the predicted values for each algorithm in one dataframe
        predDF.insert(loc=len(predDF.columns), column=libnames[counter], value= preds)
        #get performance measures (add row to resultsDF)
        resultsDF = get_performance_vals(y, preds, libnames[counter], 
                                    cv_gen, 0.95, resultsDF)
    return predDF, resultsDF

def plot_cvAUC(resultsDF, plot_title, figName, outDir, ran_Analysis):
    """
        Create plots of cvAUC along with error bars
    """
    if ran_Analysis==True:
        labels = resultsDF.index.values #array with method labels
    #hacky way to compensate for lack of index values in imported csv -- tofix
    else:
        labels = resultsDF['Unnamed: 0']

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
    plt.xlim(.5, 1) #random guessing gives .5 AUC so we need not go lower 
    fig.tight_layout()
    #style.use('ggplot') #alternative 1 -- use defaults that mimic R's ggplot
    #rstyle(ax) #alternative 2 -- gives same result as style.use('ggplot')
    plt.savefig(outDir + 'A_' + figName + '.eps', dpi=1200)
    #plt.show() 
    plt.close()
    
def plot_ROC(y, predDF, resultsDF, figName, outDir, ran_Analysis):
    """
        Plots the ROC curves using actual y and predicted probabilities
    """
    #get list of libnames, ordered by cvAUC
    if ran_Analysis==True:
        libnames = resultsDF.index.values #array with method labels
    #hacky way to compensate for lack of index values in imported csv -- tofix
    else:
        libnames = resultsDF['Unnamed: 0']
    
    plt.figure(figsize=(6.5,6.5))
    for counter, libname in enumerate(libnames):
        pred_prob = predDF[libname]
        #points to plot
        false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y, pred_prob)
        roc_auc = metrics.auc(false_positive_rate, true_positive_rate)
        #plot of points
        mycolors = get_colors()
        plt.plot(false_positive_rate, true_positive_rate, lw=2.2, 
            color=mycolors[counter], label=libname+', AUC = %0.2f'% roc_auc)
    #create legend and labels etc. and save graph
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend(loc='best')
    plt.savefig(outDir + 'C_' + figName + '.eps', dpi=1200)
    #plt.show()
    plt.close()
    

def add_info_and_print(resultsDF, include_imp_dums, screenType, pred_count, 
                        patient_sample, nobs, outName, outDir, print_results):
    resultsDF.insert(loc=len(resultsDF.columns), column='include_imp_dums',
                     value=include_imp_dums)
    resultsDF.insert(loc=len(resultsDF.columns), column='screen_method',
                     value=screenType)
    resultsDF.insert(loc=len(resultsDF.columns), column='predictor_count',
                     value=pred_count)
    resultsDF.insert(loc=len(resultsDF.columns), column='patient_sample',
                     value=patient_sample)
    resultsDF.insert(loc=len(resultsDF.columns), column='n', value=nobs)
    #sort results by decreasing cvAUC
    resultsDF.sort(columns=['cvAUC','errRate'], axis=0, ascending=False, inplace=True)
    print "\nCombined results: \n" , resultsDF
    #output results to csv
    if print_results==True:
        resultsDF.to_csv(outDir + 'R_' + outName + '.txt', sep=",")    
    #alternative: could use the predicted values coming from running SL
    #print sl.y_pred_cv.shape # numpy array of dimension n (#obs) by k (#algorithms)
    #print "\nPerformance results, RF: " , get_performance_vals(y, sl.y_pred_cv[:,0])
    return resultsDF

def get_VIM1(predictors, df, y, outName, FileNamePrefix, predictor_desc, patient_sample,
                    include_imp_dums, screenType, libs, libnames, run_VIM1, outDir):
    """
        ## Run super learner in absence of each variable separately ##
        # (could consider also dropping the corresponding imputation dummy)
        # and with all variables  -- report differences in cvAUC with full model
        #note: predictor list will include imputation dummies if that option is selected
    """
    myName = 'VIM1_' + outName
    if run_VIM1==True:
        resultsDFvim1 = pd.DataFrame() #empty df to hold performance measures
        for counter, var in enumerate(predictors):
            print "VIM1 for: " , var
            pminus1 = [p for p in predictors if p != var]
            Xvim1 = df[pminus1].astype(float).values     
            cv_gen = cv.StratifiedKFold(y, n_folds=5, shuffle=True, random_state=10)
            sl = SuperLearner(libs, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
            SL_predsvim1 = cross_val_predict_proba(sl, Xvim1, y, cv_gen)
            resultsDFvim1 = get_performance_vals(y, SL_predsvim1, var, 
                                            cv_gen, 0.95, resultsDFvim1)
        # Add results from Super Learner without dropping any variables (for comparison)
        X = df[predictors].astype(float).values
        SL_predsvim1 = cross_val_predict_proba(sl, X, y, cv_gen)
        resultsDFvim1 = get_performance_vals(y, SL_predsvim1, "Super Learner", 
                                            cv_gen, 0.95, resultsDFvim1)
        # Add column to hold varnames 
            #currently in index but that created merge problem later
        resultsDFvim1['varname'] = resultsDFvim1.index.values
        # Add columns with additional methods info, print results to text file #
        pred_count = len(predictors) - 1
        resultsDFvim1 = add_info_and_print(resultsDFvim1, include_imp_dums, screenType,
                pred_count, patient_sample, df.shape[0], myName, 
                outDir, print_results=True)
    else:
        # Read in text file if you do not wish to re-create VIM results
        resultsDFvim1 = pd.read_csv(outDir + 'R_' + myName + '.txt', sep=',')
    # Calculate VIM measure and keep only relevant variables
    # We want var importance to be 100*(cvAUC_SL_noDrops - cvAUC)
    resultsDFvim1['cvAUC'] *= 100
    cvAUC_SL_noDrops = resultsDFvim1.loc[resultsDFvim1['varname']=='Super Learner']['cvAUC'].values
    resultsDFvim1['SL_VariableDrop'] = cvAUC_SL_noDrops - resultsDFvim1['cvAUC'] 
    resultsDFvim1 = resultsDFvim1[['varname','SL_VariableDrop']]   
    return resultsDFvim1


def get_VIM2(predictors, df, y, outName, FileNamePrefix, predictor_desc, 
                include_imp_dums, screenType, patient_sample, run_VIM2, outDir):
    """
        ## Run super learner with each variable separately -- report cvAUC 
        #use SL library that only includes algorithms that work with 1 predictor
    """
    libsUni, libnamesUni = build_library( p=len(predictors), nobs=df.shape[0], 
                        screen=None, testing=False, univariate=True)
    print "univariate libnames: ", libnamesUni
    meanlib = [dummy.DummyClassifier(strategy='most_frequent')]
    myName = 'VIM2_' + outName
    if run_VIM2==True:
        resultsDFvim2 = pd.DataFrame() #empty df to hold performance measures
        for counter, var in enumerate(predictors):
            print "VIM2 for: " , var
            Xvim2 = df[[var]].astype(float).values     
            cv_gen = cv.StratifiedKFold(y, n_folds=5, shuffle=True, random_state=10)
            #cannot do prediction if not enough variation in data
            #identify binary predictors that have fewer than 10 obs in any group 
                #todo: get var list from clinic_VarsD.txt rather than manually entering here
            binary_vars = ["Dolor_Abdominal","AlterConciencia","Artralgia",
            "Ascitis","Escalofrio","ExtremidadesFrias","Tos","CianosisOP","DificultadRes",
            "Epigastralgia","Cefalea","HIPOTERM","Ictericia","Mialgia","CONGNASAL","NAUSEAS",
            "Palidez","Derrame_","Perdida_Apetito","Llenado_Capilar","PulsePressure","Rash",
            "Dolor_ocular","GargantaEnrrojecida","Sudoracion","Taquicardia","Vomito","Encias",
            "Equimosis","Nariz","Hematoma","Hematuria","Hematemesis","Hemoptisis",
            "Hipermenorrea","Melena","PetequiasEspontaneas","Purpura","Subconjuntival",
            "Vaginal","Venopuncion","Hemoconcentracion","HematocritoElev",
            "Ascitis_Signos_Radiologicos","Fluido_Para_Peri_Renal","Fluido_Perivesicular",
            "Edema_","Alveolar","Cardiomegalia","Interticial","Neumonia"]
            valcounts = df.groupby(var).size()
            if var in binary_vars and (min(valcounts) < 5 or valcounts.shape[0]==1) :
                print "Too few values for prediction - assuming no importance - run mean"
                sl = SuperLearner(meanlib, loss="nloglik", K=5, 
                                stratifyCV=True, save_pred_cv=True)    
            else:
                sl = SuperLearner(libsUni, loss="nloglik", K=5, 
                                stratifyCV=True, save_pred_cv=True)
            SL_predsvim2 = cross_val_predict_proba(sl, Xvim2, y, cv_gen)
            resultsDFvim2 = get_performance_vals(y, SL_predsvim2, var, 
                                            cv_gen, 0.95, resultsDFvim2)
        # Add column to hold varnames 
            #currently in index but that created merge problem later
        resultsDFvim2['varname'] = resultsDFvim2.index.values
        # Add columns with additional methods info, print results to text file #
        print "resultsDFvim2: " , resultsDFvim2
        pred_count = 1 #number of predictors (for output data)
        resultsDFvim2 = add_info_and_print(resultsDFvim2, include_imp_dums, screenType,
                pred_count, patient_sample, df.shape[0], myName, 
                outDir, print_results=True)
    else:
        # Read in text file if you do not wish to re-create VIM results
        resultsDFvim2 = pd.read_csv(outDir + 'R_' + myName + '.txt', sep=',')
    # Calculate VIM measure and keep only relevant variables
    resultsDFvim2['SL_Univariate'] = resultsDFvim2['cvAUC']*100
    resultsDFvim2 = resultsDFvim2[['varname','SL_Univariate']]
    return resultsDFvim2

def plot_VIMs(resultsVIM, outDir, figname):
        #print "to graph: " 
        #print resultsVIM[["varname","RF_OOB","RF_Gini","SL_Univariate","SL_VariableDrop"]]
        VIMlist = ["RF_OOB", "RF_Gini", "SL_Univariate", "SL_VariableDrop"]
        VIM_labels = ["RF - Permutation", "RF - Gini", "SL - Univariate", "SL - Variable Drop"]
        #sort by variable category and then by importance value
        resultsVIM.sort(columns=['CC_broadcat_sort','RF_OOB'], axis=0, 
            ascending=[False,True], inplace=True)
        # plot VIM results in one graph
        plt.figure(figsize=(6.9,8))
        positions = np.arange(resultsVIM.shape[0]) + .5
        #symbols and lines for for each VIM
        mymarkers = ['s','o','^','*'] 
        mymarkersizes = [5, 5 , 5, 8]
        mylstyles = ['-','-','--','-'] 
        mynoimportvals = [0, 0, 50, 0]
        #colors to indicate variable category
        colors = ["yellow","green","blue","purple"] #1st color will be ignored
        mycolors = sns.xkcd_palette(colors) #get_colors() #sns.color_palette("Set2", 10)
        clist=[mycolors[catnum] for catnum in resultsVIM['CC_broadcat_sort'].values ]  
        for counter, VIM in enumerate(VIMlist):
            #values to plot
            importances = resultsVIM[VIM]                 
            #plot of points
            plt.scatter(importances, positions, marker=mymarkers[counter],
                color=clist, label=VIM)
            #add verticle lines to indicate no importance
            plt.axvline(x=mynoimportvals[counter], linestyle = mylstyles[counter], 
                    ymin=0, ymax=resultsVIM.shape[0], color='0', linewidth=.5)
        #make left spacing large enough for labels.  Default is  .1, .9, .9, .1
        plt.subplots_adjust(left=.2, right=.9, top=.9, bottom=.1)
        #create legend and labels etc. and save graph
        plt.xlabel('Importance')
        xmin = min( [min(resultsVIM[column]) for column in VIMlist] ) - 1
        print "xmin: " , xmin
        plt.xlim(xmin,100)
        plt.ylim(0,resultsVIM.shape[0])
        plt.yticks(positions, np.array(resultsVIM["CC_name"]))
        plt.tick_params(axis='y', labelsize=6)
        #get the coloring of y-axis labels to correspond to variable cost categories
        [l.set_color(clist[i]) for i,l in enumerate(plt.gca().get_yticklabels())]   
        #remove the tick marks; they are unnecessary with the tick lines we just plotted.  
        plt.tick_params(axis="both", which="both", bottom="off", top="off",  
                        labelbottom="on", left="off", right="off", labelleft="on") 
        #create a custom legend so I can make the colors gray 
            #otherwise legend markers will be color of last variable category plotted
        lhandles = []
        for counter, VIM in enumerate(VIMlist):
            hand = mlines.Line2D([], [], fillstyle='full', color='0', linewidth=.5,
                        marker=mymarkers[counter], linestyle = mylstyles[counter],
                        markersize=mymarkersizes[counter])
            lhandles.append(hand)
        plt.legend((lhandles), (VIM_labels), prop={'size':8})
        plt.savefig(outDir + figname + '.eps', dpi=1200)
        #plt.show()
        plt.close()

def main():
    start_time_overall = time.time()

    ###########################################################################
    ###################### Choices, choices, choices ##########################

    ## Choose which parts of code to run ##
    run_MainAnalysis = False  # if false, will expect to get results from file
    plot_MainAnalysis = False  # if true, will create figures for main analysis
    run_testdata = False  # true means to get predictions for independent test set
    run_VIM = True  # if false, none of the VIM code will be run
    run_VIM1 = False # if false, will expect to obtain VIM1 results from file
    run_VIM2 = True  # if false, will expect to obtain VIM2 results from file

    ## Choose outcome variable ##
    #outcome = "is.DEN"  
    outcome = "is.DHF_DSS"

    ## Choose whether to exclude OFI patients ##
    NoOFI = True #only applies to is.DHF_DSS analyses

    ## Choose whether to exclude samples with initial DHF/DSS diagnosis ##
    NoInitialDHF = True #only applies to is.DHF_DSS analyses (should generally select True)

    ## Choose patient sample ##
    #patient_sample = "all" 
    #patient_sample = "cohort_only"
    patient_sample = "hospital_only"

    ## Choose sample to treat as independent test set (if run_testdata=True) ##
    testSample = "hospital" #options: "cohort" or "hospital
    
    ## Choose list of variables to use in prediction ##
    predictor_desc = "covarlist_all"
    #predictor_desc = "covarlist_noUltraX"
    #predictor_desc = "covarlist_CohortRestrict"
    #predictor_desc = "covarlist_genOnly"
    #predictor_desc = "covarlist_custom"  #one-off analysis

    ## Choose whether to include  ##
    include_study_dum = False #true to include is.cohort indicator 
    include_imp_dums = False #true to add imputation dummies to covariate list
    imp_dums_only = False #true to run with imputation dummies and no other variable values

    ## Choose input data (this data was prepared in R) ##
    inputData = "clin12_full_wImputedRF1.txt"
 
    ## Choose variable screening method (if any) ##
    screenType = None #current options: None, "univariate", "L1"
    screenNum = 5 #applies when screenType != None; else will be ignored

    ## Choose whether to standardize predictors (will not apply to imputation dummies)
    std_vars = True 
    
    ## Use a tiny SL library for testing purposes
    testlib = False #false if you want to run normally

    ############# you are now done with choosing your parameters ##############
    ###########################################################################


    ## Ensure consistency btwn parameter settings and develop titles for output
    if outcome=="is.DEN":
        comparison_groups = "OFI vs. DENV using " #will appear in graph title
        FileNamePrefix = "OFI.v.DEN" #use all samples; exclude PCR predictors
        NoInitialDHF = False #whether to exclude samples with initial DHF/DSS diagnosis
        NoOFI = False
    elif outcome=="is.DHF_DSS":
        if NoOFI == True:
            comparison_groups = "DF vs. DHF/DSS using " #will appear in graph title
            FileNamePrefix = "DF.v.DHFDSS"
        else:
            comparison_groups = "OFI/DF vs. DHF/DSS using " #will appear in graph title
            FileNamePrefix = "OFIDF.v.DHFDSS"
    if (NoInitialDHF==True):
        restrictions = ", DENV patients with non-severe initial Dx"
    else:
        restrictions = ""
    if patient_sample=="cohort_only" or patient_sample=="hospital_only":
        include_study_dum = False

    ## Suffix to indicate use of imputation dummies
    if imp_dums_only==True:
        FileNameSuffix = "_dumsOnly"
    elif include_imp_dums==True:
        FileNameSuffix = "_impDums"
    elif include_study_dum==True: # include is.cohort indicator
        FileNameSuffix = "_studyDum"
    else: # all covariates, but no missing indicators
        FileNameSuffix = "" #could use "_noDums" but instead I will use no suffix

    ## Preliminary list of predictors ##
    predictors_prelim = get_predictor_desc(predictor_desc+".txt", outcome, NoOFI)
    #predictors_prelim = ['Melena']  #test
    #print "Original predictor list:\n" , predictors_prelim

    ## Create pandas dataframe with data that was cleaned in R ##
    df, predictors = get_data(inputsDir, inputData, 
                    NoInitialDHF, patient_sample, NoOFI, outcome, predictors_prelim,
                    include_study_dum, include_imp_dums, imp_dums_only, standardize=std_vars)
    print "Predictors to include, pre-screening:\n" , predictors

    ## Build library of classifiers ##
    screen, pred_count = get_screen(screenType, screenNum, predictors)
    libs, libnames = build_library( p=len(predictors), nobs=df.shape[0], 
                            screen=screenType, testing=testlib)
    print "libnames: ", libnames

    ## Keep only columns in predictors list, create arrays ##
    X = df[predictors].astype(float).values
    y = df[outcome].astype(int).values #make outcome 0/1 and convert to np array
    #print "Actual outcomes: " , y[:10]
    #X, y=datasets.make_classification(n_samples=88, n_features=95) #toy data

    ## Name of text file containing performance measures (to either create or import)
    outName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + FileNameSuffix

    ## Get CV predictions and performance measures ##
    if run_MainAnalysis == True:
        cv_gen = cv.StratifiedKFold(y, n_folds=5, shuffle=True, random_state=10)
        predDF = pd.DataFrame() #this will hold predicted probs for all algorithms
        resultsDF = pd.DataFrame() #empty df to hold performance measures
        # Super Learner
        start_time_cvSL = time.time()
        sl = SuperLearner(libs, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
        SL_preds = cross_val_predict_proba(sl, X, y, cv_gen)
        predDF.insert(loc=len(predDF.columns), column='Super Learner', value=SL_preds)
        resultsDF = get_performance_vals(y, SL_preds, "Super Learner", 
                                        cv_gen, 0.95, resultsDF)
        log_statement("\ncvSL execution time: {} minutes".format(
            (time.time() - start_time_cvSL)/60. ) ) 
        # Results for each algorith in library
        predDF, resultsDF = results_for_library(X, y, cv_gen, libs, libnames, 
                            predDF, resultsDF)
        # print predicted probabilities to file (optional)
        predDF.to_csv(outDir+ 'P_' + outName + '.txt', sep=",") 
        ## Add columns with additional methods info, print results to text file ##
        resultsDF = add_info_and_print(resultsDF, include_imp_dums, screenType,
                    pred_count, patient_sample, df.shape[0], outName, 
                    outDir, print_results=True)
    else:
        print "reading from results files"
        predDF = pd.read_csv(outDir + 'P_' + outName + '.txt', sep=',') 
        resultsDF = pd.read_csv(outDir + 'R_' + outName + '.txt', sep=',') 

    if plot_MainAnalysis == True:
        ## Make bargraphs of cvAUC results ##
        plot_title = comparison_groups+predictor_desc+restrictions
        plot_cvAUC(resultsDF, plot_title="", figName=outName, 
                    outDir=outDir, ran_Analysis=run_MainAnalysis)
        ## ROC curve plots ##
        plot_ROC(y, predDF, resultsDF, outName, outDir, run_MainAnalysis)

    ## Get predictions and performance measures for test (cohort) data ##
    if run_testdata == True:
        myoName = FileNamePrefix+'_'+predictor_desc+'_'+testSample+'Test'
        dfnew = get_data(inputsDir, inputData, NoInitialDHF, 
                testSample+"_only", NoOFI, outcome, predictors, 
                include_study_dum, include_imp_dums, imp_dums_only, standardize=std_vars)
        Xnew = dfnew[predictors].astype(float).values
        ynew = dfnew[outcome].astype(int).values #make outcome 0/1 and convert to np array
        predDFnew, resultsDFnew = results_for_testset(X, y, Xnew, ynew, 
                                libs, libnames, sl_folds=5)
        # print predicted probabilities to file (optional)
        predDFnew.to_csv(outDir+ 'P_' + myName + '.txt', sep=",") 
        ## Add columns with additional methods info, print results to text file ##
        resultsDFnew = add_info_and_print(resultsDF, include_imp_dums, screenType,
                    pred_count, testSample+"Test", df.shape[0], myoName, 
                    outDir, print_results=True)
        plot_ROC(ynew, predDFnew, resultsDFnew, myoName, outDir, ran_Analysis=True)

    ## Get variable importance measures ##
        #takes about 8 min per variable when run on Nandi (~11 hr for 85 vars)
    if run_VIM==True:

        #get results from the "variable drop" importance method      
        resultsDFvim1 = get_VIM1(predictors, df, y, outName, FileNamePrefix, 
                                predictor_desc, include_imp_dums, screenType, 
                                patient_sample, libs, libnames, run_VIM1, outDir)
        #get results from the "univariate" importance method
        resultsDFvim2 = get_VIM2(predictors, df, y, outName,
                        FileNamePrefix, predictor_desc, include_imp_dums, screenType, 
                        patient_sample, run_VIM2, outDir)

        ## Read in VIMs from random forests run in R ##
            #currently just done for hospit_only, covarlist_all
        VIMout_from_R = "VIM_rf_" + FileNamePrefix + FileNameSuffix
        VIM_rf = pd.read_csv(inputsDir + VIMout_from_R + '.txt', sep='\t') 
        VIM_rf = VIM_rf.rename(columns={'variable.name.in.final.data': 'varname'})

        ## Combine VIMs into 1 dataframe ##
        prelim = pd.merge(resultsDFvim1[['varname','SL_VariableDrop']], VIM_rf, 
                            on="varname", how="inner", sort=False)
        resultsVIM = pd.merge(resultsDFvim2[['varname','SL_Univariate']], prelim, 
                            on="varname", how="inner", sort=False)

        ## Plot VIM results in one graph ##
        plot_VIMs(resultsVIM, outDir, 'VIMs_' + outName)
    
    ## Ouput execution time info ##
    log_statement("Total execution time: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()
