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
print "Version of pandas: " , pd.__version__ #should be v0.17.1
from pandas.core.categorical import Categorical

import matplotlib as mpl
mpl.use('Agg') #avoids running an X server (to avoid errors with remote runs)
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
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

inputsDir = "/srv/scratch/ccotter/intermediate_data/" #Nandi
outDir = "/users/ccotter/python_out/" #Nandi (contains small output files)
boutDir = "/srv/scratch/ccotter/py_out/" #Nandi (contains cleaned LCMS data)

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
    if (include_imp_dums==True or imp_dums_only==True):
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

def get_data(inputsDir, filename, inLCMSData, NoInitialDHF, patient_sample, 
            NoOFI, onlyLCMSpatients, noLCMSpatients, outcome, predictors_prelim, include_study_dum, 
            include_imp_dums, imp_dums_only, include_clinvars, include_LCMSvars, standardize=True):
    """
        Create pandas dataframe from text files created in R and extract_LCMS_features.py
        Standardize predictor variables (if standardize=True)
        Eliminate columns according to:
            predictors_prelim (List of covariates, excluding study and imputation indictors)
            include_study_dum (True to include is.cohort variable)
            include_imp_dums (True to include imputation dummy variables)
        Add columns according to:
            include_LCMSvars (True to add LCMS variables)
        Eliminate observations according to:
            NoInitialDHF: True means to exclude patients with initial severe dengue Dx
            patient_sample: "cohort_only","hospital_only" or "both"
            NoOFI: True means to eliminate non-dengue patients
            onlyLCMSpatients: True means to eliminate patients without LCMS data
            noLCMSpatients: True to eliminate non-LCMS patients
    """
    df_prelim = pd.read_csv(inputsDir + filename, sep='\t', 
        true_values=["True","Yes","TRUE"], false_values=["False","No","FALSE"], 
        na_values=['NaN','NA']) 

    if include_clinvars==True:

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
        nonXdf = df_prelim[['code','Study','DENV','WHO_initial_given','Cod_Nin',outcome] + addons]
        if imp_dums_only == True: #only use imputation indicators
            df = nonXdf 
        else: #include clinical info plus possible addons
            df = pd.concat([nonXdf, X], axis=1)

    else:
        #no need for standardizing and building predictor list if not including clin vars
        df = df_prelim
        all_predictors = []
    
    #restrict to non-LCMS patients if desired
    if noLCMSpatients==True:
        LCMS = pd.read_pickle(boutDir+inLCMSData)
        outerj = pd.merge(df, LCMS, on=('code','Study','Cod_Nin'), how='outer', indicator="_merge")
        df = outerj[outerj._merge=='left_only'] #keep when _merge var is left_only
    #add in LCMS data if desired 
    if include_LCMSvars==True or onlyLCMSpatients==True:
        LCMS = pd.read_pickle(boutDir+inLCMSData)
        #inner join to keep only LCMS patients
        df = pd.merge(df, LCMS, on=('code','Study','Cod_Nin'))
    #add LCMS variables if desired 
        #standardize these? They are on same scale as one another but not with clinvars
        #not a huge deal since we will generally do VIM stuff separately for clin and LCMS
    if include_LCMSvars==True:
        LCMSvars_found = []
        for var in df.columns.values:
            y = re.findall("MZ_.*", var)
            assert len(y) in (0,1)
            if len(y)==1: LCMSvars_found.append(y[0])
        all_predictors = all_predictors + LCMSvars_found
    
    #remove rows according to parameter values
    if (NoInitialDHF==True):
        df = df[df.WHO_initial_given!="DSS"] 
        df = df[df.WHO_initial_given!="DHF"] 
    if (NoOFI==True):
        df = df[df.DENV=="Positivo"] #limit to only DENV positive patients
    if (patient_sample == "hospital_only"):
        print "Run for only hospital patients"
        df = df[df.Study=="Hospital"] #limit to only hospital patients
    if (patient_sample == "cohort_only"):
        df = df[df.Study=="Cohort"] #limit to only cohort patients
    print "Number of rows in dataframe: " , df.shape[0], "\n"
    #print "Column names in dataframe: " , list(df.columns.values), "\n" 

    #useful for trouble-shooting
    #print df[["code","Study","Cod_Nin",outcome]]
    df[["code","Study","Cod_Nin",outcome]].to_csv(outDir+ 'patients_and_dx.txt', sep=",") 

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
    n_jobs = 2
    cv_grid = 5

    mean = dummy.DummyClassifier(strategy='most_frequent') #predict most frequent class

    #very simple CART (for comparison sake)
        #try max_depths of 1 thru 20
        #todo: bump min_samples_leaf up to 2 (1 is default)
            #ex: 'min_samples_leaf':[2] 
    CARTparams = {'max_depth': list(np.arange(21)[1:])} 
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

    NNparams = {'n_neighbors':[3,5,7,9,11], 'weights':['uniform','distance']}
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
        libs=[CARTune] #for testing
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

    sns.set_style("whitegrid") #necessary for getting back graph frame
    mpl.rcParams['lines.color'] = 'white'
    mpl.rcParams['text.color'] = 'white'

    #get list of libnames, ordered by cvAUC
    if ran_Analysis==True:
        libnames = resultsDF.index.values #array with method labels
    #hacky way to compensate for lack of index values in imported csv -- tofix
    else:
        libnames = resultsDF['Unnamed: 0']
    
    plt.figure(figsize=(6,6)) #(6.5,6.5) for dissertation; 
    plt.grid(b=True, which='both', axis='both',color='0.3',linestyle='-')
    for counter, libname in enumerate(libnames):
        pred_prob = predDF[libname]
        #points to plot
        false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y, pred_prob)
        #output these points to a file (optional)
        if libname in ['CART','Super Learner']:
            print counter , ": Retreiving ROC for " , libname
            d = {'fpr':pd.Series(false_positive_rate,index=thresholds),
                'tpr':pd.Series(true_positive_rate,index=thresholds)}
            df = pd.DataFrame(d)
            df.to_csv(outDir+ 'roc'+libname + '_' + figName + '.txt', sep=",") 
        roc_auc = metrics.auc(false_positive_rate, true_positive_rate)
        #plot of points
        mycolors = get_colors()
        plt.plot(false_positive_rate, true_positive_rate, lw=2.2, 
            color=mycolors[counter], label=libname+', AUC = %0.2f'% roc_auc)
    #create legend and labels etc. and save graph
    plt.xlabel('False positive rate', color='white',size=16)
    plt.ylabel('True positive rate', color='white',size=16)
    plt.xticks(color='white',size=14)
    plt.yticks(color='white',size=14)
    leg = plt.legend(loc='best', prop={'size':16})
    for text in leg.get_texts():
        plt.setp(text, color='white', size=14)
    plt.title('ROC curves for predicting servere dengue', size=16)
    #note: .png can be done with transparent background while .eps cannot
    plt.savefig(outDir + 'C_' + figName + '.png', dpi=1200, transparent=True)
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
    #sort results by decreasing cvAUC if it exists; otherwise, sort by AUC
    if patient_sample not in ['hospitalTest','cohortTest']:
        resultsDF.sort(columns=['cvAUC','errRate'], axis=0, ascending=False, inplace=True)
    else:
        resultsDF.sort(columns=['AUC','errRate'], axis=0, ascending=False, inplace=True)
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
    dumlib = [tree.DecisionTreeClassifier(min_samples_leaf=2, max_depth=2)]
    myName = 'VIM2_' + outName
    if run_VIM2==True:
        resultsDFvim2 = pd.DataFrame() #empty df to hold performance measures
        for counter, var in enumerate(predictors):
            print "VIM2 for: " , var
            Xvim2 = df[[var]].astype(float).values     
            cv_gen = cv.StratifiedKFold(y, n_folds=5, shuffle=True, random_state=10)
            #cannot do prediction if not enough variation in data
            #identify binary predictors that have fewer than Xx obs in any group 
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
            impy = re.findall(".+?_imputed", var) #length will be 1 if var has _imputed suffix
            valcounts = df.groupby(var).size()
            if (var in binary_vars or len(impy)==1):
                #LDA, QDA, trees with many branches do not make sense for binary predictors
                    #just use CART with 1 branch (and mean), provided sufficient variation
                    #otherwise, just use mean
                if (min(valcounts) < 5 or valcounts.shape[0]==1) :
                    print "Too few values for prediction - assuming no importance - run mean"
                    sl = SuperLearner(meanlib, loss="nloglik", K=5, 
                                stratifyCV=True, save_pred_cv=True)   
                else:
                    print "Use dummy library - CART and means"
                    sl = SuperLearner(dumlib+meanlib, loss="nloglik", K=5, 
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

def plot_VIMs(resultsVIM, outDir, figname, forget_VIM1):
        """
            Plot various indicators of variable importance.
            VIM measures are expected to be found in resultsVIM dataframe
            If forget_VIM1 is True, will not plot VIM1
        """

        sns.set_style("whitegrid") #necessary for getting back graph frame
        mpl.rcParams['lines.color'] = 'white'
        mpl.rcParams['text.color'] = 'white'
        

        #print "to graph: " 
        #print resultsVIM[["varname","RF_OOB","RF_Gini","SL_Univariate","SL_VariableDrop"]]
        if forget_VIM1==False:
            VIMlist = ["RF_OOB", "RF_Gini", "SL_Univariate", "SL_VariableDrop"]
            VIM_labels = ["RF - Permutation", "RF - Gini", "SL - Univariate", "SL - Variable Drop"]
        else:
            VIMlist = ["RF_OOB", "RF_Gini", "SL_Univariate"]
            VIM_labels = ["RF - Permutation", "RF - Gini", "SL - Univariate"]
        #sort by variable category and then by importance value
        resultsVIM.sort(columns=['CC_broadcat_sort','SL_Univariate'], axis=0, 
            ascending=[False,True], inplace=True)
        # plot VIM results in one graph
        plt.figure(figsize=(6.9,8)) #(6.9,8)
        plt.grid(b=True, which='both', axis='both',color='0.3',linestyle='-')
        positions = np.arange(resultsVIM.shape[0]) + .5
        #symbols and lines for for each VIM
        mymarkers = ['s','o','^','*'] 
        mymarkersizes = [5, 5 , 5, 8]
        mylstyles = ['-','-','--','-'] 
        mynoimportvals = [0, 0, 50, 0]
        #colors to indicate variable category
        #colors = ["yellow","green","blue","purple"] #1st color will be ignored
        colors = ["yellow","yellow green","cyan","pastel purple"] #JSM poster
        mycolors = sns.xkcd_palette(colors) #get_colors() #sns.color_palette("Set2", 10)
        clist=[mycolors[catnum] for catnum in resultsVIM['CC_broadcat_sort'].values ]  
        for counter, VIM in enumerate(VIMlist):
            #values to plot
            importances = resultsVIM[VIM]                 
            #plot of points
            print "importances: ", importances
            print "positions: ", positions
            plt.scatter(importances, positions, marker=mymarkers[counter],
                color=clist, label=VIM)
            #add verticle lines to indicate no importance
            plt.axvline(x=mynoimportvals[counter], linestyle = mylstyles[counter], 
                    ymin=0, ymax=resultsVIM.shape[0], color='1', linewidth=.5)
        #make left spacing large enough for labels.  Default is  .1, .9, .9, .1
        plt.subplots_adjust(left=.2, right=.9, top=.9, bottom=.1)
        #create legend and labels etc. and save graph
        plt.xlabel('Importance', color='white')
        xmin = min( [min(resultsVIM[column]) for column in VIMlist] ) - 1
        print "xmin: " , xmin
        plt.xlim(xmin,100)
        plt.xticks(color='white')
        plt.ylim(0,resultsVIM.shape[0])
        plt.yticks(positions, np.array(resultsVIM["CC_name"]))
        plt.tick_params(axis='y', labelsize=6)
        #get the coloring of y-axis labels to correspond to variable cost categories
        [l.set_color(clist[i]) for i,l in enumerate(plt.gca().get_yticklabels())]   
        #remove the tick marks; they are unnecessary with the tick lines we just plotted.  
        #plt.tick_params(axis="both", which="both", bottom="off", top="off",  
        #                labelbottom="on", left="off", right="off", labelleft="on") 
        #create a custom legend so I can make the colors gray 
            #otherwise legend markers will be color of last variable category plotted
        lhandles = []
        for counter, VIM in enumerate(VIMlist):
            hand = mlines.Line2D([], [], fillstyle='full', color='1', linewidth=.5,
                        marker=mymarkers[counter], linestyle = mylstyles[counter],
                        markersize=mymarkersizes[counter])
            lhandles.append(hand)
        plt.title('Importance of variables for distinguishing OFI from dengue')
        leg = plt.legend((lhandles), (VIM_labels), prop={'size':8})
        for text in leg.get_texts():
            plt.setp(text, color='white')
        plt.savefig(outDir + figname + '.png', dpi=1200, transparent=True)
        #plt.show()
        plt.close()

def plot_only_VIM2(resultsDFvim2):
    """
        Take a look at univariate SL results with imputation dummies
        note: R's random forests results aren't currently available for imputation dummies
            #need to fix variable naming problem when adding variable labels from cohortD 
    """
    plt.subplots_adjust(left=.3, right=.9, top=.9, bottom=.1)
    importances = resultsDFvim2['SL_Univariate']
    positions = np.arange(resultsDFvim2.shape[0]) + .5
    plt.scatter(importances, positions)
    plt.tick_params(axis='y', labelsize=6)
    plt.yticks(positions, np.array(resultsDFvim2["varname"]))
    plt.savefig(outDir + "VIM2_only" + '.eps', dpi=1200)

def subset_w_ttest(X, y, J, LCMS_indices):
    #find test stat and pval for each LCMS variable 
    tstats = []
    t_indices = []
    for i in LCMS_indices:
        s_neg = X[y==0, i]
        s_pos = X[y==1, i]
        tstat, pval = stats.ttest_ind(s_neg, s_pos, equal_var=False)
        #invoke pval criteria since high test stat is worthless if high pval
        if pval<.05:
            tstats.append(tstat)
            t_indices.append(i)
    tstats = np.array(t_indices)
    t_indices = np.array(t_indices)
    #indices in tstats of topJ LCMS features 
    topJ_tindices = np.argpartition(tstats, -J)[-J:]
    topJ_tindices = topJ_tindices[np.argsort(tstats[topJ_tindices])] #sorted
    #original indices of the J features
    topJ_indices = t_indices[topJ_tindices]
    print "values of topJ_indices: " , tstats[topJ_tindices] 
    return topJ_indices

def subset_w_topRF(X, y, J, LCMS_indices):
    #obtain array of variable importance scores
    rf = RandomForestClassifier(random_state=101)
    rf = rf.fit(X,y)
    all_VIM_rf = rf.feature_importances_
    #keep only LCMS features in ranked list
    LCMS_VIM_rf = all_VIM_rf[LCMS_indices]
    #indices in LCMS_VIM_rf of topJ LCMS features 
    topJ_nindices = np.argpartition(LCMS_VIM_rf, -J)[-J:]
    topJ_nindices = topJ_nindices[np.argsort(LCMS_VIM_rf[topJ_nindices])] #sorted
    #original indices of these chosen features
    topJ_indices = LCMS_indices[topJ_nindices]
    print "values of topJ_indices: " , all_VIM_rf[topJ_indices]
    return topJ_indices

def add_one_more_index(LCMS_indices_totry, indices_of_covars, X, y, topJ_indices):
    """Of the variables contained in LCMS_indices_totry, find one that does best
            when combined with variables in indices_of_covars.
        Return updated LCMS_indices_totry and indices_of_covars
    """         
    #cvAUC for RF done for each variable to find best var to add
    resultsDF = pd.DataFrame() #empty df to hold performance measures
    for i in LCMS_indices_totry:
        X_i = X[:, np.hstack((indices_of_covars, i)) ] 
        cv_gen = cv.StratifiedKFold(y, n_folds=3, shuffle=True, random_state=10)
        p = X_i.shape[1] #number of variables
        #print "X_i.shape " , X_i.shape, X_i
        if p < 9:
            RF_max_features = [ int(math.ceil(math.sqrt(p))) ]
        else:
            RF_max_features = [ int(math.ceil(math.sqrt(p))), int(math.ceil(p/3)) ]
        RFparams = {'n_estimators':[500],  'max_features': RF_max_features}
        #set parameters for GridSearchCV
        n_jobs = 2
        cv_grid = 3
        RFtune = grid_search.GridSearchCV(RandomForestClassifier(), RFparams,
                score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)
        preds = cross_val_predict_proba(RFtune, X_i, y, cv_gen)[:, 1]
        #the row index of resultsDF will equal i, the index that we care about 
        resultsDF = get_performance_vals(y, preds, str(i), cv_gen, .95, resultsDF)
    #find index of best LCMS variable
    row_max = np.argmax(np.array(resultsDF['cvAUC'])) #row with highest cvAUC
    index = int(resultsDF.index.tolist()[row_max]) #row index of this row
    #print "resultsDF" , resultsDF
    #print "row_max: ", row_max
    print "index: " , index
    #update which indices to try and which to include as covars for next round
    LCMS_indices_totry = LCMS_indices_totry[np.where(LCMS_indices_totry!=index)]
    indices_of_covars = np.hstack((indices_of_covars, [index]))
    topJ_indices.append(index)
    return LCMS_indices_totry, indices_of_covars, topJ_indices

def get_topJ(X, y, predictors, J, subsetAlg):
    """Return the positions (in predictor list) of top LCMS predictors as numpy array
       subsetAlg specified desired method ('topRF' or 'ttest' etc.)
    """
    #find index numbers of LC-MS features
    LCMS_indices = []
    nonLCMS_indices = []
    for counter, name in enumerate(predictors):
        f = re.findall("MZ_.*", name)
        assert len(f) in (0,1)
        if len(f)==1: LCMS_indices.append(counter)
        elif len(f)==0: nonLCMS_indices.append(counter)
    LCMS_indices = np.array(LCMS_indices)
    nonLCMS_indices = np.array(nonLCMS_indices)

    if subsetAlg=="ttest":
        topJ_indices = subset_w_ttest(X, y, J, LCMS_indices) 

    elif subsetAlg=="topRF":
        topJ_indices = subset_w_topRF(X, y, J, LCMS_indices)

    elif subsetAlg=="greedyRF":
        topJ_indices = []
        LCMS_indices_totry = LCMS_indices
        indices_of_covars = nonLCMS_indices
        #add another index J times to get top J indices
        for j in xrange(J):
            LCMS_indices_totry, indices_of_covars, topJ_indices = add_one_more_index(
                LCMS_indices_totry, indices_of_covars, X, y, topJ_indices)
        topJ_indices = np.array(topJ_indices)    

    print "topJ_indices: " , topJ_indices
    
    return topJ_indices

def TopJ_and_clin(X, predictors, topJ_indices, j):
    """Input:
            "X" is a numpy array with columns as listed in "predictors" list
            Indices of top J LCMS predictors are in "topJ_indices" list
       Output:
            (1) a new X numpy array containing only the clin data and 
               'j' LCMS columns, which are the top j of the topJ_indices
            (2) a new predictors list to correspond to the new X numpy array
    """

    #find indices of non-LCMS variables
    nonLCMS_indices = []
    for counter, name in enumerate(predictors):
        y = re.findall("MZ_.*", name)
        assert len(y) in (0,1)
        if len(y)==0: nonLCMS_indices.append(counter)
    nonLCMS_indices = np.array(nonLCMS_indices)
    print "nonLCMS_indices: " , nonLCMS_indices

    #add in the indices of top j LCMS variables
    if j>0:
        print "topJ_indices: " , topJ_indices
        desired_indices = np.concatenate((nonLCMS_indices, topJ_indices[-j:]), axis=0)
    else:
        desired_indices = nonLCMS_indices
    print "desired_indices: " , desired_indices

    #now rebuild X and predictors list based on indice list
    X_new = X[ : , desired_indices]
    predictors_new = [predictors[i] for i in desired_indices]
    print "predictors_new: " , predictors_new

    return X_new, predictors_new

def predict_for_LCMSpatients_allClin(inputsDir, inClinData, inLCMSData,
            NoInitialDHF, NoOFI, outcome, predictors_prelim, 
            include_study_dum, include_imp_dums, imp_dums_only, std_vars, testlib):
    """Predictions for LCMS patients using clin data fit to all non-LCMS hospit patients
            no need to keep track of predictor indices and whatnot
    """

    patient_sample = "hospital_only"

    #fit SL with all non-LCMS observations
    onlyLCMSpatients = False
    noLCMSpatients = True
    include_clinvars = True
    include_LCMSvars = False
    df_noLCMS, predictors_clin = get_data(inputsDir, inClinData, inLCMSData,
            NoInitialDHF, patient_sample, NoOFI, onlyLCMSpatients, noLCMSpatients,
            outcome, predictors_prelim, include_study_dum, include_imp_dums, imp_dums_only, 
            include_clinvars, include_LCMSvars, std_vars)
    libs, libnames = build_library( p=len(predictors_clin), nobs=df_noLCMS.shape[0],
                     screen=None, testing=testlib)
    sl_noLCMS = SuperLearner(libs, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
    X_noLCMS = df_noLCMS[predictors_clin].astype(float).values
    print "type and X_noLCMS.shape: " , type(X_noLCMS), X_noLCMS.shape
    y_noLCMS = df_noLCMS[outcome].astype(int).values 
    print "type and y_noLCMS.shape: " , type(y_noLCMS), y_noLCMS.shape
    clinSL = sl_noLCMS.fit(X_noLCMS, y_noLCMS)

    #get predicted probs for LCMS data using fitted model from above
       #no need to keep track of predictor indices and whatnot - using only clin data 
    onlyLCMSpatients = True
    noLCMSpatients = False
    include_clinvars = True
    include_LCMSvars = False
    df_onlyLCMS, predictors_clin = get_data(inputsDir, inClinData, inLCMSData,
            NoInitialDHF, patient_sample, NoOFI, onlyLCMSpatients, noLCMSpatients,
            outcome, predictors_prelim, include_study_dum, include_imp_dums, imp_dums_only, 
            include_clinvars, include_LCMSvars, standardize=std_vars)
    X_onlyLCMS = df_onlyLCMS[predictors_clin].astype(float).values
    print "type and X_onlyLCMS.shape: " , type(X_onlyLCMS), X_onlyLCMS.shape
    y_onlyLCMS = df_onlyLCMS[outcome].astype(int).values 
    print "type and y_onlyLCMS.shape: " , type(y_onlyLCMS), y_onlyLCMS.shape
    SL_preds_clin = clinSL.predict_proba(X_onlyLCMS)

    return y_onlyLCMS, SL_preds_clin

def bestSubset_predictions(inputsDir, inClinData, inLCMSData,
            NoInitialDHF, NoOFI, outcome, predictors_prelim, 
            include_study_dum, include_imp_dums, imp_dums_only, std_vars, testlib,
            predsSub, resultsSub, SL_preds_clin, subsetMethod, J, subsetAlg):
    """predictions for LCMS patients using topJ LCMS features plus 
        (a) clinical vars -- subsetMethod = 1
        (b) SL_preds_clin -- subsetMethod = 2 
        (c) both clin vars and SL_preds_clin -- subsetMethod = 3
       select topJ LCMS features using covars that are also used in prediction
    """
    #obtain data for LCMS obs with predictor set (=LCMS data plus various clin vars)
    onlyLCMSpatients = True
    noLCMSpatients = False
    include_LCMSvars = True
    if subsetMethod==1 or subsetMethod==3:
        include_clinvars = True
    elif subsetMethod==2:
        include_clinvars = False
    df_allLCMS, predictors_all = get_data(inputsDir, inClinData, inLCMSData,
            NoInitialDHF, "both", NoOFI, onlyLCMSpatients, noLCMSpatients,
            outcome, predictors_prelim, include_study_dum, include_imp_dums, imp_dums_only, 
            include_clinvars, include_LCMSvars, std_vars)
    #add SL_preds_clin as a covariate if desired
    if subsetMethod==2 or subsetMethod==3:
        df_allLCMS["SL_preds_clin"] = SL_preds_clin
        predictors_all = predictors_all + ["SL_preds_clin"]
    #turn pandas dataframe into data arrays 
    X_allLCMS = df_allLCMS[predictors_all].astype(float).values
    print "type and X_allLCMS.shape: " , type(X_allLCMS), X_allLCMS.shape
    y_allLCMS = df_allLCMS[outcome].astype(int).values 
    print "type and y_allLCMS.shape: " , type(y_allLCMS), y_allLCMS.shape

    #obtain topJ LCMS predictors using data created above
    topJ_indices = get_topJ(X_allLCMS, y_allLCMS, predictors_all, J, subsetAlg)

    #obtain predictions using same data fed to get_topJ, but with only topJ LCMS vars
    cv_gen = cv.StratifiedKFold(y_allLCMS, n_folds=3, shuffle=True, random_state=10)
    LCMS_var_count = sum( [len(re.findall("MZ_", var)) for var in predictors_all] )
    print LCMS_var_count, " LCMS features in data."
    for j in range(J) + [J] + [LCMS_var_count]:
    #for j in [J]: #testing
        print "******* Running for top " , j , "LCMS vars, method " , subsetMethod , " *******"
        if j < LCMS_var_count:
            X_reduced, predictors_reduced = TopJ_and_clin(
                X_allLCMS, predictors_all, topJ_indices, j)
        else:
            #data and predictor list for no LCMS restrictions
            X_reduced = X_allLCMS
            predictors_reduced = predictors_all
        print "X_reduced.shape: " , X_reduced.shape
        libs, libnames = build_library( p=len(predictors_reduced), nobs=df_allLCMS.shape[0],
                            screen=None, testing=testlib)
        sl_topJ = SuperLearner(libs, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
        SL_preds_topJ = cross_val_predict_proba(sl_topJ, X_reduced, y_allLCMS, cv_gen)
        print "SL_preds_topJ: " , SL_preds_topJ.shape, SL_preds_topJ
        colname = 'M' + str(subsetMethod) + '_LCMS_' + str(j)
        predsSub.insert(loc=len(predsSub.columns), column=colname, value=SL_preds_topJ)
        resultsSub = get_performance_vals(y_allLCMS,SL_preds_topJ,colname,cv_gen, 0.95,resultsSub)

    return predsSub, resultsSub

def bestSubset_analysis(inputsDir, inClinData, inLCMSData,
                NoInitialDHF, NoOFI, outcome, predictors_prelim, 
                include_study_dum, include_imp_dums, imp_dums_only, 
                std_vars, testlib, J, subsetAlg):
    """Output predicted probabilities and performance results csv
        for various algorithms with up to J LCMS features chosen
    """

    predsSub = pd.DataFrame() #this will hold predicted probs for all runs here
    resultsSub = pd.DataFrame() #empty df to hold performance measures
    
    #predictions for LCMS patients using clin data fit to all non-LCMS hospit patients
        #no need to keep track of predictor indices and whatnot
    y_onlyLCMS, SL_preds_clin = predict_for_LCMSpatients_allClin(
                inputsDir, inClinData, inLCMSData,
                NoInitialDHF, NoOFI, outcome, predictors_prelim, 
                include_study_dum, include_imp_dums, imp_dums_only, std_vars, testlib)

    #predictions for LCMS patients using topJ LCMS features
    for method in [1, 2, 3]: # [1,2,3] to do for all 3 subset methods
        predsSub, resultsSub = bestSubset_predictions(inputsDir, inClinData, inLCMSData,
                NoInitialDHF, NoOFI, outcome, predictors_prelim, 
                include_study_dum, include_imp_dums, imp_dums_only, std_vars, testlib,
                predsSub, resultsSub, SL_preds_clin, method, J, subsetAlg)
    print "predsSub: " , predsSub
    print "resultsSub: " , resultsSub
    
    # print predicted probabilities to file (optional)
    predsSub.to_csv(outDir+ 'P_' + outcome + '_BestSubset_' + subsetAlg + '_NPserum.txt', sep=",") 
    ## Add columns with additional methods info, print results to text file ##
    resultsSub.to_csv(outDir+ 'R_' + outcome + '_BestSubset_' + subsetAlg + '_NPserum.txt', sep=",")

    return 

def run_ttests(X,y,permute):
    tstats = []
    pvals = []
    for i in xrange(X.shape[1]):
        if permute == True:
            #permutation will be different for each feature.
            y = np.random.permutation(y)
        s_neg = X[y==0, i]
        s_pos = X[y==1, i]
        tstat, pval = stats.ttest_ind(s_neg, s_pos, equal_var=True)
        tstats.append(tstat)
        pvals.append(pval)
    tstats_ordered = np.sort(np.array(tstats))
    pvals = np.array(pvals)
    return tstats_ordered, pvals

def get_permute_stats(X, y, K):
    """get K permuted t-stats and p-vals for each LCMS feature
    """
    for i in xrange(K):
        tstats_ordered_perm, pvals_perm = run_ttests(X, y, permute=True)
        if i == 0:
            tstats_allK = tstats_ordered_perm 
            pvals_allK = pvals_perm       
        else:
            tstats_allK = np.vstack((tstats_allK, tstats_ordered_perm))
            pvals_allK = np.hstack((pvals_allK, pvals_perm))
    #ordered stat i is the average ordered stat i across the K values 
    stats_ordered_permute = np.mean(tstats_allK, axis=0)  
    return pvals_allK, tstats_allK, stats_ordered_permute

def find_delta_for_SAM(stats_ordered_permute, tstats_ordered, tstats_allK,
                     FDR, K, scaleMe):
    """find the delta that gives desired FDR
    """
    #have to make adjustments for kink when using scaleMe
    if scaleMe==True:
        #take out the region corresponding to imposed zeros
        cut_start = np.searchsorted(tstats_ordered, 0)
        cut_stop = np.searchsorted(tstats_ordered, 0, side='right')
        margin = 300
        stats_ordered_permute = np.hstack((stats_ordered_permute[:cut_start-margin],
                                           stats_ordered_permute[cut_stop+margin:]))
        tstats_ordered = np.hstack((tstats_ordered[:cut_start-margin],
                                           tstats_ordered[cut_stop+margin:]))
    tstat_diff = stats_ordered_permute - tstats_ordered
    print "stats_ordered_permute: " , stats_ordered_permute
    print "tstats_ordered: " , tstats_ordered
    print "tstat_diff: " , tstat_diff
    deltas_to_try = np.arange(.05, 1.7, .05)
    for i, delta in enumerate(deltas_to_try):
        print "Delta: " , delta
        called_bool = np.abs(tstat_diff) > delta
        called = tstats_ordered[called_bool==True, ]
        smallest_pos = min(called[called>0,])
        print "smallest_pos: " , smallest_pos
        largest_neg = max(called[called<0,])
        print "largest_neg: " , largest_neg
        #number of features called significant 
        called_count = np.sum(called_bool)
        print "called_count: " , called_count
        #number of falsely significant genes under the null
        false_call1 = stats_ordered_permute[stats_ordered_permute < largest_neg,]
        false_call2 = stats_ordered_permute[stats_ordered_permute > smallest_pos,]
        false_call_count = (len(false_call1) + len(false_call2)) 
        print "false_call_count: " , false_call_count
        #FDR is number of false positives under null divided by number called significant
        FDR_est = false_call_count / float(called_count)
        print "FDR_est: " , FDR_est
        if FDR_est < FDR:
            print "Found a winner"
            return delta
        elif i == len(deltas_to_try)-1:
            print "No winner found"
            return .4  

def scale_stats(stats_ordered_permute, tstats_ordered, X, nfeatures_orig):
    """Create faux data to compensate for pre-screening of features by lab
    """

    n_to_add = nfeatures_orig - len(tstats_ordered)
    df = X.shape[0] - 2 #degrees of freedom for t dist
  
    #create additional tstats (draw from t distribution)
        #want to have nfeatures_orig - nfeatures_postScreen
    for i in xrange(10):
        tstats_ordered_perm = np.sort(np.random.standard_t(df, n_to_add))
        if i == 0:
            tstats_allK = tstats_ordered_perm       
        else:
            tstats_allK = np.vstack((tstats_allK, tstats_ordered_perm))
    #ordered stat i is the average ordered stat i across the K values 
    stats_to_add = np.mean(tstats_allK, axis=0) 

    #combine with the tstats from permutated data
    new_stats_permute = np.sort(np.hstack((stats_ordered_permute, stats_to_add)))

    #as a "worst case", assume t-stat was 0 for features that were dropped
    ind = np.searchsorted(tstats_ordered, 0) #location to insert
    new_tstats_ordered = np.insert(tstats_ordered, ind, np.zeros(n_to_add))
    
    return new_tstats_ordered, new_stats_permute


def make_hist_of_tstats(tstats_ordered, stats_ordered_permute, 
                        outDir, outcome, inLCMSData, nfeatures_orig, scaleMe=False):
    """Histogram of t-stats using input:
        tstats_ordered are tstats based on data
        stats_ordered_permute are tstats based on permutation of data
        nfeatures_orig is number of features in data before pre-screening
    """ 
    tstat_hist = plt.figure(figsize=(4,4))
    bins = np.linspace(-5, 5, 40)
    #plot actual t-stats
    #re-norm t-stats so sum of area is (nfeatures_postScreen/nfeatures_orig) rather than 1
    if inLCMSData=="MassHuntNP" and scaleMe==True:
        scaleFactor = len(tstats_ordered) / float(nfeatures_orig)
        print "Scaling tstat histogram by factor of " , scaleFactor
        hist, bins = np.histogram(tstats_ordered, bins, normed=True)
        plt.bar(bins[:-1], hist*scaleFactor, width=bins[1]-bins[0], color='aqua')
    else:
        plt.hist(tstats_ordered, bins, normed=True, color='aqua', alpha=.5)
    #for null distribution, plot all t-stats from permutation (=K*744)
    plt.hist(stats_ordered_permute, bins, normed=True, histtype='step', 
            color='purple', lw=1)
    plt.axvline(x=0, ymin=0, ymax=10,  color='gold', linestyle='-', lw=1)
    #for another null distribution, plot the normal dist (~t-dist)
    x = np.linspace(-6,6,100)
    sigma = math.sqrt(np.var(tstats_ordered)/len(tstats_ordered))
    print "sigma is " , sigma
    plt.plot(x, mlab.normpdf(x, 0, 1), color="grey", lw=1)
    #set graph parameters
    plt.xlim(-6, 6)
    if scaleMe==True:
        #need close-up if we have scaled down values
        plt.ylim(0,.05)
        suffix = "_scaled"
    else:
        suffix = ""
    plt.xlabel("t-statistics")
    #make left spacing large enough for labels.  Default is  .1, .9, .9, .1
    plt.subplots_adjust(left=.12, right=.9, top=.9, bottom=.15)
    plt.savefig(outDir + outcome + '_LCMS_tstat_hist_' + inLCMSData + suffix + '.eps', dpi=1200) 
    plt.close() 

def make_SAM_plot(stats_ordered_permute, tstats_ordered, best_delta, 
                outDir, outcome, inLCMSData, scaleMe=False):
    """plot of expected ordered t-stats vs. actual ordered t-stats (SAM)
    """
    SAM_scatter = plt.figure(figsize=(6.7,6.7))
    #based on permuted data
    plt.scatter(stats_ordered_permute, tstats_ordered, s=3, alpha=.5)
    #based on t distribution
    #stats_ordered_tdist = np.sort(np.random.standard_t(df=86, size=len(tstats_ordered)))
    #plt.scatter(stats_ordered_tdist, tstats_ordered, s=3, alpha=.5, color='g')
    #45 degree line
    xvals = np.linspace(-4, 4, 10000)
    plt.plot(xvals, xvals, 'g-', lw=.5, alpha=.5)
    #display SAM threshold band
    plt.plot(xvals, xvals - best_delta, 'r-', lw=.5, alpha=.5)
    plt.plot(xvals, xvals + best_delta, 'r-', lw=.5, alpha=.5)
    plt.xlabel("Expected order statistics")
    plt.ylabel("t-statistic")
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    if scaleMe==True:
        suffix = "_scaled"
    else:
        suffix = ""
    #make left spacing large enough for labels.  Default is  .1, .9, .9, .1
    plt.subplots_adjust(left=.15, right=.9, top=.9, bottom=.15)
    plt.savefig(outDir + outcome + '_tstat_SAM_plot_' + inLCMSData + suffix + '.eps', dpi=1200) 
    plt.close()

def plot_BH(pvals, outDir, FDR, outcome, inLCMSData, nfeatures_orig, scaleMe=False):
    """Plot ordered pvals with Benjamini-Hochberg line

    """
    pval_scatter = plt.figure(figsize=(4,4))
    plt.scatter(np.array(range(len(pvals))), np.sort(pvals), 
        alpha=.5, s=3, edgecolors='face')
    #BH threshold is FDR*j/#features
    xvals = np.linspace(0,len(pvals),1000)
    if inLCMSData=="MassHuntNP" and scaleMe==True:
        plt.plot(xvals, xvals*FDR/float(nfeatures_orig), 'r--', lw=1)
    else:
        plt.plot(xvals, xvals*FDR/float(len(pvals)), 'r--', lw=1)
    #intersection of BH threshold and ordered pvals
    #plt.axvline(x=318, ymin=-10, ymax=10,  color='gray', lw=1)
    #uniform distribution
    plt.plot(xvals, xvals/len(pvals), 'g-', lw=1, alpha=.5)
    plt.xlabel("Molecular features ordered by p-value")
    plt.ylabel("p-value")
    if scaleMe==True:
        if outcome=="is.DHF_DSS":
            plt.xlim(-1,50) 
            plt.ylim(-.000005,.0005) 
        else:
            plt.xlim(-5,200)
            plt.ylim(-.0001,.002)
        suffix = "_scaled"
    else:
        plt.xlim(-10,760)
        plt.ylim(-.05,1.05)
        suffix = ""
    #make left spacing large enough for labels.  Default is  .1, .9, .9, .1
    plt.subplots_adjust(left=.15, right=.9, top=.9, bottom=.15)
    plt.savefig(outDir + outcome + '_pval_BH_plot_' + inLCMSData + suffix + '.eps', dpi=1200) 
    plt.close() 

def reality_check(pvals_allK, outDir):
    """verify pvals from permutation are uniformly distributed btwn 0 and 1
    """
    pval_null = plt.figure()
    plt.scatter(np.array(range(len(pvals_allK))), np.sort(pvals_allK),s=1)
    xvals = np.linspace(0,len(pvals_allK),1000)
    plt.plot(xvals, xvals/len(pvals_allK), 'g-', lw=1, alpha=.5)
    plt.savefig(outDir + 'pval_null_plot.eps', dpi=1200)
    plt.close()


def FDR_analysis(FDR, K, scaleMe, nfeatures_orig, inputsDir, inClinData, inLCMSData,
            NoInitialDHF, patient_sample, NoOFI, onlyLCMSpatients, noLCMSpatients,
            outcome, predictors_prelim, 
            include_study_dum, include_imp_dums, imp_dums_only, std_vars, outDir):   
    """Run several FDR related analyses and produce corresponding graphs

    """
    #obtain LCMS data with diagnostic outcome (but not other covariates)
    include_clinvars = False 
    include_LCMSvars = True
    df, predictors = get_data(inputsDir, inClinData, inLCMSData,
                NoInitialDHF, patient_sample, NoOFI, onlyLCMSpatients, noLCMSpatients,
                outcome, predictors_prelim, 
                include_study_dum, include_imp_dums, imp_dums_only, 
                include_clinvars, include_LCMSvars, standardize=std_vars)
    X_orig = df[predictors].astype(float).values
    X = preprocessing.scale(X_orig) #standardized values
    print "type and X.shape: " , type(X), X.shape
    y = df[outcome].astype(int).values 
    print "type and y.shape: " , type(y), y.shape        

    #obtain actual t-stats and p-vals from data for each LCMS feature
    tstats_ordered, pvals = run_ttests(X, y, permute=False)

    #plot of ordered pvals with BH threshold line
    plot_BH(pvals, outDir, FDR, outcome, inLCMSData, nfeatures_orig, scaleMe)

    #get K permuted t-stats and p-vals for each LCMS feature
    pvals_allK, tstats_allK, stats_ordered_permute = get_permute_stats(X, y, K)

    #histogram of t-stats
    make_hist_of_tstats(tstats_ordered, stats_ordered_permute, 
            outDir, outcome, inLCMSData, nfeatures_orig, scaleMe)

    #find the delta that gives desired FDR for SAM method and make SAM plot
    if inLCMSData=="MassHuntNP" and scaleMe==True:
        #create modified tstat arrays to account for # of originally extracted features
        tstats_ordered, stats_ordered_permute = scale_stats(
                stats_ordered_permute, tstats_ordered, X, nfeatures_orig)
    best_delta = find_delta_for_SAM(stats_ordered_permute, tstats_ordered, 
                                tstats_allK, FDR, K, scaleMe)    
    #plot of expected ordered t-stats vs. actual ordered t-stats (SAM)
    make_SAM_plot(stats_ordered_permute, tstats_ordered, best_delta, 
            outDir, outcome, inLCMSData, scaleMe)

    #curiousity -- do we get uniform dist of pvals for permutation? yes :)
    #reality_check(pvals_allK, outDir)


def parse_arguments():
    """Parse arguments provided at the command line level
    """        

    import argparse
    parser = argparse.ArgumentParser(description='Predict dengue.')
    
    parser.add_argument('--run_MainAnalysis', action='store_true', default=False)
    parser.add_argument('--run_SL', action='store_true', default=False)
    parser.add_argument('--plot_MainAnalysis', action='store_true', default=False)
    parser.add_argument('--run_testdata', action='store_true', default=False)
    parser.add_argument('--run_BestSubset', action='store_true', default=False)
    parser.add_argument('--run_VIM', action='store_true', default=False)
    parser.add_argument('--run_VIM1', action='store_true', default=False)
    parser.add_argument('--forget_VIM1', action='store_true', default=False)
    parser.add_argument('--run_VIM2', action='store_true', default=False)

    parser.add_argument('--outcome', dest='outcome')

    parser.add_argument('--NoOFI', action='store_true', default=False)
    parser.add_argument('--NoInitialDHF', action='store_true', default=False)
    
    parser.add_argument('--patient_sample', dest='patient_sample',
                        default = "hospital_only")
    parser.add_argument('--testSample', dest='testSample')

    parser.add_argument('--include_clinvars', action='store_true', default=False)
    parser.add_argument('--include_LCMSvars', action='store_true', default=False)
    parser.add_argument('--onlyLCMSpatients', action='store_true', default=False)

    parser.add_argument('--predictor_desc', dest='predictor_desc',
                        default = "covarlist_all")

    parser.add_argument('--include_study_dum', action='store_true', default=False)
    parser.add_argument('--include_imp_dums', action='store_true', default=False)
    parser.add_argument('--imp_dums_only', action='store_true', default=False)

    parser.add_argument('--inClinData', dest='inClinData', 
                        default = "clin12_full_wImputedRF1.txt")
    parser.add_argument('--inLCMSData', dest='inLCMSData', default='')

    parser.add_argument('--testlib', action='store_true', default=False)

    args = parser.parse_args()

    return (args.run_MainAnalysis, args.run_SL, args.plot_MainAnalysis, args.run_testdata, 
            args.run_BestSubset, args.run_VIM, args.run_VIM1, args.run_VIM2,
            args.outcome,
            args.NoOFI, args.NoInitialDHF, 
            args.patient_sample, args.testSample, 
            args.include_clinvars, args.include_clinvars, args.include_LCMSvars, 
            args.onlyLCMSpatients, args.predictor_desc,
            args.include_study_dum, args.include_imp_dums, args.imp_dums_only,
            args.inClinData, args.inLCMSData, args.testlib)

def main():
    start_time_overall = time.time()

    ###########################################################################
    ###################### Choices, choices, choices ##########################
    
    #set to true if running "run_prediction_in_python.py loops thru parameter lists
    params_from_commandline = False 

    if params_from_commandline==True:
        ## Parse parameters provided at command-line level
        (run_MainAnalysis, run_SL, plot_MainAnalysis, run_testdata, 
         run_BestSubset, run_VIM, run_VIM1, run_VIM2,
         outcome, NoOFI, NoInitialDHF, patient_sample, testSample, 
         include_clinvars, include_clinvars, include_LCMSvars, onlyLCMSpatients,
         predictor_desc,
         include_study_dum, include_imp_dums, imp_dums_only, 
         inClinData, inLCMSData, testlib
         ) = parse_arguments()

    else:
        ## Choose which parts of code to run ##
        run_MainAnalysis = False  # if false, will expect to get results from file
        run_SL = False # if false, can still run algorithms in library, but not SL
        plot_MainAnalysis = False  # if true, will create figures for main analysis
        run_testdata = False  # true means to get predictions for independent test set
        run_BestSubset = True # true means to find best subset of LCMS features and run SL with them
        run_FDR = False #true means to do false discovery rate analysis
        #determine which variable importance code to run
        run_VIM = False  # if false, none of the VIM code will be run
        forget_VIM1 = False # if true, will do VIM analysis and graphs without VIM1 output
        run_VIM1 = False # if false, will expect to obtain VIM1 (multivariate) results from file
        run_VIM2 = False  # if false, will expect to obtain VIM2 (univariate) results from file
        only_VIM2 = False # true if just want to plot VIM2 (and not other VIMs)

        ## Choose outcome variable ##
        #outcome = "is.DEN"  
        outcome = "is.DHF_DSS"

        ## Choose whether to exclude OFI patients ##
        NoOFI = False #only applies to is.DHF_DSS analyses

        ## Choose whether to exclude samples with initial DHF/DSS diagnosis ##
        NoInitialDHF = False #only applies to is.DHF_DSS analyses (should generally select True)

        ## Choose patient sample - only applies when not using LCMS data ##
        #patient_sample = "all" 
        #patient_sample = "cohort_only"
        patient_sample = "hospital_only"

        ## Choose sample to treat as independent test set (if run_testdata=True) ##
        testSample = "hospital" #options: "cohort" or "hospital"

        ## Choose whether to include clinical variables and/or LCMS features ##
        include_clinvars = True #true to include clinical variables
        include_LCMSvars = True #true to include LCMS features

        ## Choose whether to restrict to only LCMS patients ##
        onlyLCMSpatients = True #true to include only patients with LCMS data
        #eventually, instead have trainObs="all"/"LCMSonly" and testObs="all"/"LCMSonly"        

        ## Choose list of clinical variables to use in prediction 
            #applicable if include_clinvars==True (else will be ignored)
        predictor_desc = "covarlist_all"
        #predictor_desc = "covarlist_noUltraX"
        #predictor_desc = "covarlist_CohortRestrict"
        #predictor_desc = "covarlist_genOnly"
        #predictor_desc = "covarlist_custom"  #one-off analysis

        ## Choose whether to include various indicators ##
        include_study_dum = False #true to include is.cohort indicator 
        include_imp_dums = False #true to add imputation dummies to covariate list
        imp_dums_only = False #true to run with imputation dummies and no other variable values

        ## Choose input clinical data ##
        inClinData = "clin12_full_wImputedRF1.txt" #data prepared in R

        ## Choose input LCMS data ##
        #inLCMSData = "NPbins50x50" #prepared in extract_LCMS_features.py
        #inLCMSData = "RPbins50x50" #prepared in extract_LCMS_features.py
        inLCMSData = "MassHuntNP" #prepared in prepare_MassHunter_data.py
        #inLCMSData = "SalivaMH" #prepared in prepare_MassHunter_data.py 
        #inLCMSData = "UrineMH" #prepared in prepare_MassHunter_data.py 

        ## Use a tiny SL library for testing purposes
        testlib = False #false if you want to run normally
     
    ## Choose variable screening method (if any) ##
    screenType = None #current options: None, "univariate", "L1"
    screenNum = 5 #applies when screenType != None; else will be ignored
    #choose True to include clinical covariates when selecting best LCMS features
        #only applies when screenType!=None, include_clinvars=T and include_LCMSvars=T

    ## Choose whether to standardize predictors (will not apply to imputation dummies)
    std_vars = True 

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
    #if we are restricting to LCMS patients we will not want to restrict on hospital/cohort
    if include_LCMSvars==True or onlyLCMSpatients==True:
        patient_sample = "all"
        print "patient sample set to all"
    ## Suffix to indicate use of imputation dummies
    if imp_dums_only==True:
        FileNameSuffix = "_dumsOnly"
    elif include_imp_dums==True:
        FileNameSuffix = "_impDums"
    elif include_study_dum==True: # include is.cohort indicator
        FileNameSuffix = "_studyDum"
    else: # all covariates, but no missing indicators
        FileNameSuffix = "" #could use "_noDums" but instead I will use no suffix
    #will be reset as needed in code below (used for LCMS best subset analysis)
    noLCMSpatients = False

    ## Check for illegal parameter setting combos
    try:
        assert include_clinvars==True or include_LCMSvars == True
    except AssertionError:
        print "Parameter error: At least one of include_clinvars or include_LCMSvars must be True"
        exit(1)

    ###########################################################################

    ## Preliminary list of predictors ##
    predictors_prelim = get_predictor_desc(predictor_desc+".txt", outcome, NoOFI)
    if True == False:
        #predictors_prelim = ['Melena']  #test
        #print "Original predictor list:\n" , predictors_prelim

        ## Create pandas dataframe with data that was cleaned in R ##
        df, predictors = get_data(inputsDir, inClinData, inLCMSData,
                        NoInitialDHF, patient_sample, NoOFI, onlyLCMSpatients, noLCMSpatients,
                        outcome, predictors_prelim, 
                        include_study_dum, include_imp_dums, imp_dums_only, 
                        include_clinvars, include_LCMSvars, standardize=std_vars)
        #print "Predictors to include, pre-screening:\n" , predictors

        ## Build library of classifiers ##
        screen, pred_count = get_screen(screenType, screenNum, predictors)
        libs, libnames = build_library( p=len(predictors), nobs=df.shape[0], 
                                screen=screenType, testing=testlib)
        print "libnames: ", libnames

        ## Keep only columns in predictors list, create arrays ##
        X = df[predictors].astype(float).values
        print "type and X.shape: " , type(X), X.shape
        y = df[outcome].astype(int).values #make outcome 0/1 and convert to np array
        print "type and y.shape: " , type(y), y.shape
        #print "Actual outcomes: " , y[:10]
        #X, y=datasets.make_classification(n_samples=88, n_features=95) #toy data

    #sl = SuperLearner(libs, loss="nloglik", K=2, stratifyCV=True, save_pred_cv=True)
    #TestSL = sl.fit(X, y)

    ## Establish name of text file containing performance measures (to either create or import)
        #used for main analysis and VIM analysis
    if onlyLCMSpatients == True: 
        #indicate restriction to only LCMS patients
            #could potentially only use clinical vars with these patients (keep predictor_desc)
        outName = FileNamePrefix + '_' + predictor_desc + '_'+inLCMSData+'patients' + FileNameSuffix
    if include_LCMSvars == True:
        #note: use of LCMS features implies restriction to LCMS patients
        if include_clinvars == False:
            outName = FileNamePrefix + '_covarlist_' + inLCMSData + FileNameSuffix
        else:
            outName = FileNamePrefix + '_' + predictor_desc + '_' + inLCMSData + FileNameSuffix
    if onlyLCMSpatients == False and include_LCMSvars == False:
        #case in which LCMS data is not used at all
        outName = FileNamePrefix + '_' + predictor_desc + '_' + patient_sample + FileNameSuffix
    print "outName: " , outName 

    ## Analysis of LCMS feature significance (FDR analysis) ##
    if run_FDR == True:
        #set parameter values
        FDR = .20
        K = 1000 #number of permutations to run (about 30 min for K=10,000) 
        scaleMe = True #true for calcs to do worst case scenario based on original #features
        nfeatures_orig = 15930 #15930 was reported as # features before filtering (MassHuntNP)
        #run analyses and create FDR related graphs
        FDR_analysis(FDR, K, scaleMe, nfeatures_orig, inputsDir, inClinData, inLCMSData,
                    NoInitialDHF, patient_sample, NoOFI, onlyLCMSpatients, noLCMSpatients,
                    outcome, predictors_prelim,
                    include_study_dum, include_imp_dums, imp_dums_only, std_vars, 
                    outDir)            

    ## Analysis with chosen subsets of LCMS features ##
    if run_BestSubset == True:
        #Output csv files with predicted probabilities and performance results
        J = 5 #select up to J top LCMS features
        subsetAlg = "greedyRF" #options: "topRF","ttest", "greedyRF"
        bestSubset_analysis(inputsDir, inClinData, inLCMSData,
                    NoInitialDHF, NoOFI, outcome, predictors_prelim, 
                    include_study_dum, include_imp_dums, imp_dums_only, 
                    std_vars, testlib, J, subsetAlg)

    ## Get CV predictions and performance measures ##
    if run_MainAnalysis == True:
        cv_gen = cv.StratifiedKFold(y, n_folds=5, shuffle=True, random_state=10)
        predDF = pd.DataFrame() #this will hold predicted probs for all algorithms
        resultsDF = pd.DataFrame() #empty df to hold performance measures
        # Super Learner
        start_time_cvSL = time.time()
        if run_SL == True:
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
        if True==True:
            alg_list = ['CART', 'Centroids', 'LDA+shrinkage', 'Gradient Boost', 
                        'AdaBoost', 'Random Forests', 'N.Neighbor', 'Logit-L1', 
                        'Logit-L2', 'SVM-L2','Super Learner']
            rDF = pd.read_csv(outDir + 'R_' + outName + '.txt', sep=",")
            resultsDF = rDF[rDF['Unnamed: 0'].isin(alg_list)]
            predDF = pd.read_csv(outDir + 'P_' + outName + '.txt', sep=',')
        else:
            predDF = pd.read_csv(outDir + 'P_' + outName + '.txt', sep=',') 
            resultsDF = pd.read_csv(outDir + 'R_' + outName + '.txt', sep=',') 

    if plot_MainAnalysis == True:
        ## Make bargraphs of cvAUC results ##
        #plot_title = comparison_groups+predictor_desc+restrictions
        #plot_cvAUC(resultsDF, plot_title="", figName=outName, 
        #            outDir=outDir, ran_Analysis=run_MainAnalysis)
        ## ROC curve plots ##
        plot_ROC(y, predDF, resultsDF, outName, outDir, run_MainAnalysis)

    ## Get predictions and performance measures for test (cohort) data ##
    if run_testdata == True:
        myoName = FileNamePrefix+'_'+predictor_desc+'_'+testSample+'Test'
        dfnew, predictors = get_data(inputsDir, inClinData, inLCMSData, NoInitialDHF, 
                testSample+"_only", NoOFI, onlyLCMSpatients, noLCMSpatients, outcome, predictors, 
                include_study_dum, include_imp_dums, imp_dums_only, 
                include_clinvars, include_LCMSvars, standardize=std_vars)
        Xnew = dfnew[predictors].astype(float).values
        ynew = dfnew[outcome].astype(int).values #make outcome 0/1 and convert to np array
        predDFnew, resultsDFnew = results_for_testset(X, y, Xnew, ynew, 
                                libs, libnames, sl_folds=5)
        # print predicted probabilities to file (optional)
        predDFnew.to_csv(outDir+ 'P_' + myoName + '.txt', sep=",") 
        ## Add columns with additional methods info, print results to text file ##
        resultsDFnew = add_info_and_print(resultsDFnew, include_imp_dums, screenType,
                    pred_count, testSample+"Test", df.shape[0], myoName, 
                    outDir, print_results=True)
        plot_ROC(ynew, predDFnew, resultsDFnew, myoName, outDir, ran_Analysis=True)

    ## Get variable importance measures ##
        #takes about 8 min per variable when run on Nandi (~11 hr for 85 vars)
    if run_VIM==True:

        #get results from the "variable drop" importance method      
        if forget_VIM1==False:
            resultsDFvim1 = get_VIM1(predictors, df, y, outName, FileNamePrefix, 
                                    predictor_desc, include_imp_dums, screenType, 
                                    patient_sample, libs, libnames, run_VIM1, outDir)
        #get results from the "univariate" importance method
        resultsDFvim2 = get_VIM2(predictors, df, y, outName,
                        FileNamePrefix, predictor_desc, include_imp_dums, screenType, 
                        patient_sample, run_VIM2, outDir)
        if only_VIM2 == True:
            plot_only_VIM2(resultsDFvim2)
        else:
            ## Read in VIMs from random forests run in R ##
                #currently just done for hospit_only, covarlist_all
            VIMout_from_R = "VIM_rf_" + FileNamePrefix + FileNameSuffix
            VIM_rf = pd.read_csv(inputsDir + VIMout_from_R + '.txt', sep='\t')
            print VIMout_from_R
            print "VIM_rf: " , VIM_rf 
            VIM_rf = VIM_rf.rename(columns={'variable.name.in.final.data': 'varname'})

            ## Combine VIMs into 1 dataframe ##
            if forget_VIM1==False:
                prelim = pd.merge(resultsDFvim1[['varname','SL_VariableDrop']], VIM_rf, 
                                    on="varname", how="inner", sort=False)
            else:
                prelim = VIM_rf
            print "resultsDFvim2: " ,  resultsDFvim2
            resultsVIM = pd.merge(resultsDFvim2[['varname','SL_Univariate']], prelim, 
                                on="varname", how="inner", sort=False)
            print "resultsVIM: " , resultsVIM

            ## Plot VIM results in one graph ##
            plot_VIMs(resultsVIM, outDir, 'VIMs_' + outName, forget_VIM1)
    
    ## Ouput execution time info ##
    log_statement("Total execution time: {} minutes".format(
        (time.time() - start_time_overall)/60. ) ) 

if __name__ == '__main__':
    main()
