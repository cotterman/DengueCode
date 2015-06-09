
from sklearn.cross_validation import *

#import warnings
from itertools import chain, combinations
from math import ceil, floor, factorial
#import numbers
#import times
from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp

from sklearn.base import is_classifier, clone
from sklearn.utils import indexable, check_random_state, safe_indexing
from sklearn.utils.validation import (_is_arraylike, _num_samples,
                               check_array, column_or_1d)
from sklearn.utils.multiclass import type_of_target
from sklearn.externals.joblib import Parallel, delayed, logger
from sklearn.externals.six import with_metaclass
from sklearn.externals.six.moves import zip
from sklearn.metrics.scorer import check_scoring
from sklearn.utils.fixes import bincount

from sklearn.cross_validation import _check_cv, _safe_split, _check_is_partition


def cross_val_predict_proba(estimator, X, y=None, cv=None, n_jobs=1,
                      verbose=0, fit_params=None, pre_dispatch='2*n_jobs'):
    """Generate cross-validated estimates for each input data point
    Parameters
    ----------
    estimator : estimator object implementing 'fit' and 'predict'
        The object to use to fit the data.
    X : array-like
        The data to fit. Can be, for example a list, or an array at least 2d.
    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.
    cv : cross-validation generator or int, optional, default: None
        A cross-validation generator to use. If int, determines
        the number of folds in StratifiedKFold if y is binary
        or multiclass and estimator is a classifier, or the number
        of folds in KFold otherwise. If None, it is equivalent to cv=3.
        This generator must include all elements in the test set exactly once.
        Otherwise, a ValueError is raised.
    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.
    verbose : integer, optional
        The verbosity level.
    fit_params : dict, optional
        Parameters to pass to the fit method of the estimator.
    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:
            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs
            - An int, giving the exact number of total jobs that are
              spawned
            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'
    Returns
    -------
    probs : ndarray
        This is the result of calling 'predict_proba'
    """
    X, y = indexable(X, y)

    cv = _check_cv(cv, X, y, classifier=is_classifier(estimator))
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    probs_blocks = parallel(delayed(_fit_and_predict_proba)(clone(estimator), X, y,
                                                      train, test, verbose,
                                                      fit_params)
                            for train, test in cv)
    p = np.concatenate([p for p, _ in probs_blocks])
    locs = np.concatenate([loc for _, loc in probs_blocks])
    if not _check_is_partition(locs, X.shape[0]):
        raise ValueError('cross_val_predict_proba only works for partitions')
    probs = p.copy()
    probs[locs] = p
    return probs


def _fit_and_predict_proba(estimator, X, y, train, test, verbose, fit_params):
    """Fit estimator and predict values for a given dataset split.
    Parameters
    ----------
    estimator : estimator object implementing 'fit' and 'predict'
        The object to use to fit the data.
    X : array-like of shape at least 2D
        The data to fit.
    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.
    train : array-like, shape (n_train_samples,)
        Indices of training samples.
    test : array-like, shape (n_test_samples,)
        Indices of test samples.
    verbose : integer
        The verbosity level.
    fit_params : dict or None
        Parameters that will be passed to ``estimator.fit``.
    Returns
    -------
    probs : sequence
        Result of calling 'estimator.predict_proba'
    test : array-like
        This is the value of the test parameter
    """
    # Adjust length of sample weights
    fit_params = fit_params if fit_params is not None else {}
    fit_params = dict([(k, _index_param_value(X, v, train))
                      for k, v in fit_params.items()])

    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, _ = _safe_split(estimator, X, y, test, train)

    if y_train is None:
        estimator.fit(X_train, **fit_params)
    else:
        estimator.fit(X_train, y_train, **fit_params)
    probs = estimator.predict_proba(X_test)
    return probs, test
