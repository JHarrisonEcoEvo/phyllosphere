#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 16:07:46 2021

@author: jharrison
"""

import numpy as np
import pandas as pd

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor

X = pd.read_csv("./processedData/imputed_scaled_ITS_metadata.csv")
 
############
# QC step ##
############

#Simulate data that should be predictive of diversity
#to ensure our model and diagnostic code is correctly specified.
#Code works well.
# testFeature = np.random.normal(loc = 2, scale = 1, size = len(X['shannonsISD']))

# X['testFeature'] = X['shannonsISD']*testFeature

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

Xtest_array = np.array(strat_test_set.loc[:, strat_test_set.columns != 'shannonsISD'])
Ytest_array = np.array(strat_test_set['shannonsISD']).reshape(-1,1)

# define the model
model = RandomForestRegressor(n_estimators=200,
                              criterion = 'mse',
                              min_samples_split = 2,
                              )
model.fit(X=X_array[:,10:], y=Y_array.ravel())

yhat = model.predict(Xtest_array[:,10:])

import matplotlib.pyplot as plt
plt.scatter(Ytest_array, yhat, )

plt.xlabel("Expected values")
plt.ylabel("Observed values");
t = plt.text(0.15, 0.5, "R2 = 0.6%",weight='bold')
t.set_bbox(dict(facecolor='red', alpha=0.8, edgecolor='red'))

plt.savefig("./visuals/RandomForest_shannons.pdf")

from scipy.stats import pearsonr
from sklearn.metrics import r2_score

r2_score(Ytest_array,yhat)

pearsonr(np.squeeze(Ytest_array), np.squeeze(yhat))

# evaluate the model
cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
n_scores = cross_val_score(model, X_array[:,10:], Y_array,
                           scoring='neg_mean_absolute_error', 
                           cv=cv, 
                           n_jobs=-1, 
                           error_score='raise')
# report performance
print('MAE: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))

##############
#SHAP values
############

import shap 
shap.initjs()

# Create Tree Explainer object that can calculate shap values
explainer = shap.TreeExplainer(model)

shap_values = explainer.shap_values(X = X_array[:,10:])
list_of_labels = list(strat_train_set.loc[:, strat_train_set.columns != 'shannonsISD'].columns)
feature_names=list_of_labels[10:]

shap.summary_plot(shap_values, 
                  X_array[:,10:],
                  feature_names=feature_names)

#Extract SHAP values and then remove some features and try rerunning model
shapValues = pd.DataFrame(shap_values, columns=feature_names, 
                          index = strat_train_set.loc[:,'samplename'])

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
duds = shapValues.sum(axis = 0)[abs(shapValues.sum(axis = 0)) < 0.2]

[duds.index]

import re
l = duds.index.tolist()
r = re.compile(r'.*\s+.*')
newlist = list(filter(r.match, l))
#print(newlist)

#Remove taxa
for i in newlist:
    l.remove(i)

#Now we have a list, l, that has features to be removed. Lets remove them and rerun the model 
###########################################
# Running model again after removing useless features #
###########################################

X = pd.read_csv("./processedData/imputed_scaled_ITS_metadata.csv")

#Remove features that we have deemed not helpful
for i in l:
    X.drop(i, inplace=True, axis=1)

#Do train/test split#
split = StratifiedShuffleSplit(n_splits=1, test_size=(0.2), random_state=(666))
for train_index, test_index in split.split(X, X["EN"]):
    strat_train_set = X.loc[train_index]
    strat_test_set = X.loc[test_index]

strat_test_set.shape  
strat_train_set.shape

X_array = np.array(strat_train_set.loc[:, strat_train_set.columns != 'shannonsISD'])
Y_array = np.array(strat_train_set['shannonsISD']).reshape(-1,1)

Xtest_array = np.array(strat_test_set.loc[:, strat_test_set.columns != 'shannonsISD'])
Ytest_array = np.array(strat_test_set['shannonsISD']).reshape(-1,1)

# define the model
model = RandomForestRegressor(n_estimators=200,
                              criterion = 'mse',
                              min_samples_split = 2,
                              )
model.fit(X=X_array[:,10:], y=Y_array.ravel())

yhat = model.predict(Xtest_array[:,10:])

r2 = r2_score(Ytest_array,yhat)

pearsonr(np.squeeze(Ytest_array), np.squeeze(yhat))

import matplotlib.pyplot as plt
plt.scatter(Ytest_array, yhat, )
 
plt.xlabel("Expected values")
plt.ylabel("Observed values");


t = plt.text(0.15, 0.5,
             "R2 = {}%".format(round(r2*100)),
             weight='bold')
t.set_bbox(dict(facecolor='red', alpha=0.8, edgecolor='red'))



