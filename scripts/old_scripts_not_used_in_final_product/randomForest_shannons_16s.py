#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:07:46 2021

@author: jharrison
"""

import numpy as np
import pandas as pd

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor

#feature engineering was done in the neural net program
X = pd.read_csv("../processedData/imputed_scaled_16S_metadata.csv")

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

###############################
# hyperparameter optimization #
###############################
#Code modified from: https://towardsdatascience.com/hyperparameter-tuning-the-random-forest-in-python-using-scikit-learn-28d2aa77dd74

from sklearn.model_selection import RandomizedSearchCV

n_estimators = [int(x) for x in np.linspace(start = 200, stop = 500, num = 10)]

# Number of features to consider at every split
max_features = ['auto', 'sqrt']

# Maximum number of levels in tree
max_depth = [int(x) for x in np.linspace(4, 50, num = 2)]
max_depth.append(None)

# Minimum number of samples required to split a node
min_samples_split = [2, 5, 10]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 4]
# Method of selecting samples for training each tree
bootstrap = [True, False]
# Create the random grid
random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}

from pprint import pprint
pprint(random_grid)

{'bootstrap': [True, False],
               'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
 'max_features': ['auto', 'sqrt'],
 'min_samples_leaf': [1, 2, 4],
 'min_samples_split': [2, 5, 10],
 'n_estimators': [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]}

##############
# model time #
##############

# define the model
model = RandomForestRegressor(n_estimators=200,
                              criterion = 'mse',
                              min_samples_split = 2,
                              random_state = 666
                              )
 
# Random search of hyper parameters, using 3 fold cross validation, 
# search across 100 different combinations, and use all available cores
model_random = RandomizedSearchCV(estimator = model, 
                               param_distributions = random_grid, 
                               n_iter = 100, 
                               cv = 3, 
                               verbose=2, 
                               random_state=666, 
                               n_jobs = -1) #-1 means use all processers


model_random.fit(X=X_array[:,10:], y=Y_array.ravel())

#see what we got
model_random.best_params_

yhat = model_random.predict(Xtest_array[:,10:])
from sklearn.metrics import r2_score

r2 = r2_score(Ytest_array,yhat)
r2

#save the model '
import joblib

#save model
#joblib.dump(model_random, './models/randomforest_Shannons_16s')

#To load use this syntax
#model_random = joblib.load("PATH")

#Convert the RandomizedSearchCV object to a randomForestRegressor object
tuned_model = RandomForestRegressor(**model_random.best_params_)

#sanity check
# yhat = model_random.predict(Xtest_array[:,10:])
tuned_model.fit(X=X_array[:,10:], y=Y_array.ravel())
yhat2 = tuned_model.predict(Xtest_array[:,10:])
r2 = r2_score(Ytest_array,yhat2)
r2

# import matplotlib.pyplot as plt
# plt.scatter(yhat, yhat2)

#################
# QC relic code #
#################

#Compare fit of OG model and tuned model
# def evaluate(model, test_features, test_labels):
#     predictions = model.predict(test_features)
#     errors = abs(predictions - test_labels)
#     mape = 100 * np.mean(errors / test_labels)
#     accuracy = 100 - mape
#     print('Model Performance')
#     print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
#     print('Accuracy = {:0.2f}%.'.format(accuracy))
    
#     return accuracy


# model.fit(X=X_array[:,10:], y=Y_array.ravel())
# base_accuracy = evaluate(model, Xtest_array[:,10:], Ytest_array)

# base_accuracy_tuned = evaluate(model_random, Xtest_array[:,10:], Ytest_array)


#More model plotting and evaluation
# import matplotlib.pyplot as plt
# plt.scatter(Ytest_array, yhat, )

# plt.xlabel("Expected values")
# plt.ylabel("Observed values");
# plt.text(0.1, 0.9, "R2 = {}%".format(round(r2*100)),weight='bold')

# plt.savefig("./visuals/RandomForest_shannons16s.pdf")

# from scipy.stats import pearsonr


# pearsonr(np.squeeze(Ytest_array), np.squeeze(yhat))

# # evaluate the model
# cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
# n_scores = cross_val_score(model, X_array[:,10:], Y_array,
#                            scoring='neg_mean_absolute_error', 
#                            cv=cv, 
#                            n_jobs=-1, 
#                            error_score='raise')
# # report performance
# print('MAE: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))

##############
#SHAP values
############

import shap 
shap.initjs()

# Create Tree Explainer object that can calculate shap values
explainer = shap.TreeExplainer(tuned_model)

shap_values = explainer.shap_values(X = X_array[:,10:])
list_of_labels = list(strat_train_set.loc[:, strat_train_set.columns != 'shannonsISD'].columns)
feature_names=list_of_labels[10:]

# shap.summary_plot(shap_values, 
#                   X_array[:,10:],
#                   feature_names=feature_names)

# #Extract SHAP values and then remove some features and try rerunning model
shapValues = pd.DataFrame(shap_values, columns=feature_names, 
                          index = strat_train_set.loc[:,'samplename'])

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

#write SHAP values to disk. 
shapValues.to_csv(path_or_buf=("./results/shap_values_randomforest_Shannons_16s"))

###############################################
#Feature removal following SHAP value creation
###############################################

#Figure out unimportant features
duds = shapValues.sum(axis = 0)[abs(shapValues.sum(axis = 0)) < 0.2]

#Here I am searching for taxa in this list of useless features, bc
#I WANT to keep all taxa

import re
l = duds.index.tolist()
r = re.compile(r'.*\s+.*')
newlist = list(filter(r.match, l))
#print(newlist)

#Remove taxa from list of features to remove. Thus keeping taxa.
for i in newlist:
    l.remove(i)

#Now we have a list, l, that has features to be removed. Lets remove them and rerun the model 
###########################################
# Running model again after removing useless features #
###########################################

#Note that the reason to do this is because if one has a bunch of shit features
#they will get picked during feature bagging while training the model, thus
#diluting the effect of the better features. 


X = pd.read_csv("./processedData/imputed_scaled_16S_metadata.csv")

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

# refit the model. NOTE that I am not bothering to retune the model as I think
#This goes into diminishing returns

tuned_model.fit(X=X_array[:,10:], y=Y_array.ravel())

yhat = tuned_model.predict(Xtest_array[:,10:])

r2_reduced = r2_score(Ytest_array,yhat)

r2 - r2_reduced

