
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
from scipy.interpolate import spline, interp1d
import scipy.ndimage as ndi
import statsmodels.api as sm

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

sys.path.append('/home/ccotter/temp_Dengue_code/SuPyLearner/supylearner')
import core #this is the main SuPyLearner code
from core import SuperLearner, cv_superlearner
from cross_val_utils import cross_val_predict_proba


np.random.seed(100)

inputsDir = "/home/ccotter/dengue_data_and_results_local/intermediate_data/" #home PC
outDir = "/home/ccotter/dengue_data_and_results_local/python_out/" #home PC
#inputsDir = "/home/nboley/Desktop/dengue_data_and_results_local/intermediate_data/" #home PC diff login
#inputsDir = "/srv/scratch/ccotter/intermediate_data/" #mitra
#outDir = "/srv/scratch/ccotter/python_out/" #mitra


###############################################################################
#################### Experiment with reading in data  #########################
###############################################################################

df = pd.read_csv(inputsDir + "clin12_full_wImputedRF1.txt", sep='\t', 
        true_values=["True","Yes"], false_values=["False","No"], na_values=['NaN','NA']) 
    #df["is.DEN"] = df["is.DEN"].astype(int) #ML funcs give predicted probs only for int outcomes
df[["is.female"]] = (df[["Sexo"]]=="female").astype(int) 
df.groupby('Mialgia').size() #gives counts for each value of Mialgia
df_sub = df[['Eosi','Dolor_ocular','DaysSick','Vaginal','is.DEN','code']]
df_sub.head()

#create numpy array and standardize data (subtract mean, divide by std)
    #keep outcome variable as is
X_prelim = df_sub.astype(float).values
X = preprocessing.scale(X_prelim)  
X[:5,] #view

#this does the same thing while working within the dataframe object
df_sub.dtypes
std_df = df_sub.drop('code').apply(lambda x: (x - x.mean()) / x.std() )
std_df.head() #same result as with X, above
#put back outcome variable pre-scaling
std_df[['is.DEN']] = df_sub[['is.DEN']]
std_df[['code']] = df_sub[['code']]
std_df.head() #perfecto


f = open(inputsDir + "covarlist_all.txt", 'r')
predictors = f.read().splitlines()



###############################################################################
#################### Experiment with predicted probabilities ##################
###############################################################################

X, y=datasets.make_classification(n_samples=88, n_features=95) #ge
print "Actual outcomes: " , y[:10]

p = 95
n_jobs = 5
cv_grid = 5
RFparams = {'n_estimators':[1000],  
            'max_features':[ int(math.ceil(math.sqrt(p))), int(math.ceil(p/3)) ]} 
#est = grid_search.GridSearchCV(RandomForestClassifier(), RFparams,
#        score_func=metrics.roc_auc_score, n_jobs = n_jobs, cv = cv_grid)

est=AdaBoostClassifier() 

est.fit(X,y)
cv_gen = cv.StratifiedKFold(y, n_folds=2, shuffle=True, random_state=10)
preds = cross_val_predict_proba(est, X, y, cv_gen)[:, 1]
print "Predicted proba: " , preds[:10]


###############################################################################
#################### Experiment with variable importance ######################
###############################################################################

resultsDFvim1 = pd.read_csv(outDir + 'R_VIM1_OFI.v.DEN_covarlist_all_hospital_only.txt', sep=',')
resultsDFvim2 = pd.read_csv(outDir + 'R_VIM2_OFI.v.DEN_covarlist_all_hospital_only.txt', sep=',')


# random forests - out-of-bag permutation method

# random forests - gini-coefficient improvement 

# univariate

# L1 penalty logistic

# forward step-wise

# MARS

# NSC

## built-in VIM

#this does not seem to be quite what I want....
forest = ExtraTreesClassifier(n_estimators=250, random_state=0)
forest.fit(X,y)
# this importance measure 
importances = forest.feature_importances_ #numpy array with one entry per variable
# forest.estimators_ contains a collection of n_estimators fitted trees
    # so here we take the std of feature_importances across all trees in fitted forest
std = np.std([tree.feature_importances_ for tree in forest.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]

# Print the feature ranking
print("Feature ranking:")
for f in range(10):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))







###############################################################################
#################### Experiment with ROC curve smoothing ######################
###############################################################################

def get_colors():
    # These are the "Tableau 20" colors as RGB.  
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
      
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
    for i in range(len(tableau20)):  
        r, g, b = tableau20[i]  
        tableau20[i] = (r / 255., g / 255., b / 255.)

    return tableau20

mycolors = get_colors()  
  

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#points to plot
false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(y, preds)
roc_auc = metrics.auc(false_positive_rate, true_positive_rate)
#plot of points
plt.plot(false_positive_rate, true_positive_rate, lw=2.5, 
        color=tableau20[3], label='AUC = %0.2f'% roc_auc)
#plt.show()

#smoothed plot
#fpr_smooth = np.linspace(
#    false_positive_rate.min(),false_positive_rate.max(),300)
#this spline method did not turn out
#tpr_smooth = spline(false_positive_rate,true_positive_rate,
#                    fpr_smooth)
#plt.plot(fpr_smooth, tpr_smooth, color=tableau20[5]) 
#plt.show()

#this moving average is not bad
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    #add something to clean up ends
    y_smooth[-box_pts:]=1 #need to do better than this
    return y_smooth
plt.plot(false_positive_rate, smooth(true_positive_rate,8),
        color=tableau20[5], lw=2.5, label='smooth')
plt.legend(loc='best')

#to improve further, might want to check out the R package pROC
    #If method="binormal", a linear model is fitted to the quantiles
    # of the sensitivities and specificities.
    # Smoothed sensitivities and specificities are then generated from  
    # this model on n points. 

#this interpolation method worked but is not better than raw data
    #still need to figure out how to get a best fit smooth curve
#f1 = interp1d(false_positive_rate, true_positive_rate)
#plt.plot(fpr_smooth, f1(fpr_smooth), color=tableau20[1], lw=2.5) 

#different interpolation
#yhat = savitzky_golay(y, 51, 3)
#plt.plot(false_positive_rate, yhat, color=tableau20[4], lw=2.5)

#lowess looks terrible, too
#lowess = sm.nonparametric.lowess(
#    false_positive_rate, true_positive_rate, frac=0.1)
#plt.plot(lowess[:, 0], lowess[:, 1], color=tableau20[4], lw=2.5)

#polynomial fit (did not work
#pf = np.polyfit(false_positive_rate, true_positive_rate, deg=2)
#pf_y = pf[0]+pf[1]*false_positive_rate+pf[2]*false_positive_rate**2
#plt.plot(false_positive_rate, pf_y, color=tableau20[2], lw=2.5) 
plt.show()




