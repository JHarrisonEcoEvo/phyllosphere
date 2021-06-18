#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 16:07:46 2021

@author: jharrison
"""

import sklearn
import numpy as np
import pandas as pd

from sklearn.datasets import make_regression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor

X = pd.read_csv("./processedData/imputed_scaled_ITS_metadata.csv")

num_features = []
categorical_features = []

for i in X:
    if X[i].dtype.kind in 'iufc':
        num_features.extend([i])
    else:
        categorical_features.extend([i])
        
#####################
#Do train/test split#
#####################
#Code from Geron. 

#Handy code for stratified test/train splitting with a random number 
#generator used for reproducibility

#Stratifying by compartment. May need to revisit this if some other feature
#is really important.
from sklearn.model_selection import StratifiedShuffleSplit

split = StratifiedShuffleSplit(n_splits=1, test_size=(0.2), random_state=(666))
for train_index, test_index in split.split(X, X["EN"]):
    strat_train_set = X.loc[train_index]
    strat_test_set = X.loc[test_index]

strat_test_set.shape  
strat_train_set.shape

X_array = np.array(strat_train_set.loc[:, strat_train_set.columns != 'shannonsISD'])
Y_array = np.array(strat_train_set['shannonsISD']).reshape(-1,1)


# define the model
model = RandomForestRegressor()
# evaluate the model
cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
n_scores = cross_val_score(model, X_array[:,10:], Y_array,
                           scoring='neg_mean_absolute_error', 
                           cv=cv, 
                           n_jobs=-1, 
                           error_score='raise')
# report performance
print('MAE: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))
