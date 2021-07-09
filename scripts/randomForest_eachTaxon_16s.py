#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:07:46 2021

@author: jharrison
"""

#Trains a random forest with hyper parameter optimization for each line
#in a file that provides a list of taxa to model. 

#NOTE: if you get this warning: WARN: OMP_NUM_THREADS=None =>
#... If you are using openblas if you are using openblas set OMP_NUM_THREADS=1 or risk subprocess calls hanging indefinitely
#Then set this env. variable: export OMP_NUM_THREADS=1

#Also to get this to work on Teton I had to load miniconda and make a new
#conda environment where I could install modules without problems. 
#I also had to add the conda forge channel to find certain modules.
#For reference, commands were: 
#module load  python/3.8.7 #perhaps not stricly neccessary to make env. but YOU DO want to make sure you are using
#the version of python you want before installing stuff or you will get a bunch of stupid errors
#because of pythonic compatibility woes
#module load miniconda3
#conda create --name py38 python=3.8
#conda config --add channels conda-forge
#conda activate py38
#conda install pandas -y
#conda install numpy -y #Do NOT load from the module list on TEton as it will cause you to switch back to python 2.7 and f things up.
#conda install sklearn -y
# conda install hyperopt -y
# conda install git pip #possibly needed to install hpsklearn from github without having issues due to user conflics between conda and pip
#  git clone https://github.com/hyperopt/hyperopt-sklearn.git
# #cd into cloned rep and then:
# pip3 install -e . --user
# conda install joblib -y
# conda install shap -y
# conda install re

import numpy as np
import pandas as pd
#Nonsense to force print entirety of a data frame
# pd.options.display.max_columns = None
# pd.options.display.max_rows = None
from hyperopt import tpe

from hpsklearn import HyperoptEstimator,random_forest_regression
# from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import r2_score #,make_scorer, mean_squared_error
from sklearn.model_selection import StratifiedShuffleSplit

# from sklearn.model_selection import cross_val_score, KFold, StratifiedShuffleSplit, RandomizedSearchCV
# from sklearn.utils import check_random_state
import joblib
import shap
import re
from csv import writer
import sys 

#This wrapping in an if statement specifying __main__ was added
#because of an error. I think hypsklearn was importing portions of this script
#or some other script and so when run from the CLI this script was not seen as
#the main script and was defined but not run. See this page for more:
#https://www.geeksforgeeks.org/what-does-the-if-__name__-__main__-do/  
if __name__ == '__main__':
    #Pull integer from command line that will correspond to a line in 
    #file with otu names to make models for. 
    #Doing this way so that I can choose which OTUs to model more easily
    #See the script, "choosing_taxa_to_model_with_rf.py"
    
    focal_integer = int(sys.argv[1]) - 1 #because of Pythonic indexing, subtract one from value
    
    eval_runs = sys.argv[2]
    
    possibles = pd.read_csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv")
    focal_taxon = possibles['taxa_16s'][focal_integer]
    
    print(focal_taxon)
    
    #feature engineering was done in the neural net program
    X = pd.read_csv("./processedData/imputed_scaled_16S_metadata.csv")
    
    #Bring in taxon proportions, these were wrangled in R after being generated by CNVRG
    taxa = pd.read_csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv")
    #taxa = pd.read_csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv")
    
    #Subset taxa to just sample and the taxon taken from the command line
    taxa = taxa[['sample',focal_taxon]]
    
    # Remove diversity from features
    #X.columns.values
    X = X.loc[:,X.columns != 'shannonsISD']
    #Remember that the first 9 fields are not of interest. 
    
    #Need to ensure the samples are in the correct order
    #First a bit of wrangling to get a matchable sample name field
    #We will merge X and taxa instead of reording one of them
    X['compartment'] = X['EN'].map({1: 'EN', 0: 'EP'})
    #X['compartment'].value_counts()
    X['sample'] = X['plant.x'].str.cat(X['compartment'], sep='_')
    
    X_taxa = pd.merge(X, taxa, on='sample')
    
    #####################
    #Do train/test split#
    #####################
    
    #Figure out which indices correspond to the samples that had
    #the most of the taxon. These top 100 samples will be used for a 
    #stratified sampling, since they will include the samples that actually
    #had the particular microbe. 
    
    #Determine samples with most of microbe and then make a one hot encoded 
    #feature to stratify by
    N = 100
    res = sorted(range(len(X_taxa[focal_taxon])), key = lambda sub: X_taxa[focal_taxon][sub])[-N:]
    
    X_taxa['focal_taxon_onehot'] = 0
    X_taxa.loc[res, 'focal_taxon_onehot'] = 1
    
    #Code from Geron. 
    
    #Handy code for stratified test/train splitting with a random number 
    #generator used for reproducibility
    #To do a nested stratification of compartment and the presence of focal taxon,
    #I am pasting the two together intoa a new feature. As others have said online
    #this is a brutal hack. 
    
    X_taxa['brutal_hack'] = X_taxa["EN"] + X_taxa['focal_taxon_onehot']
    
    split = StratifiedShuffleSplit(n_splits=1, test_size=(0.3), random_state=(666))
    for train_index, test_index in split.split(X_taxa, X_taxa['brutal_hack']):
        strat_train_set = X_taxa.loc[train_index]
        strat_test_set = X_taxa.loc[test_index]
    
    #define response/predictor split
    #First drop various book keeping columns
    testsamples = strat_test_set['samplename']
    trainsamples = strat_train_set['samplename']
    
    strat_test_set = strat_test_set.drop(['focal_taxon_onehot', 'brutal_hack', 'sample',
     'region_site',
     'compartment',
     'Unnamed: 0',
     'X',
     'taxon.x',
     'region_site_plant',
     'plant.x',
     'forward_barcode',
     'reverse_barcode',
     'locus',
     'samplename'], axis=1)
    
    strat_train_set = strat_train_set.drop(['focal_taxon_onehot', 'brutal_hack', 'sample',
     'region_site',
     'Unnamed: 0',
     'compartment',
     'X',
     'taxon.x',
     'region_site_plant',
     'plant.x',
     'forward_barcode',
     'reverse_barcode',
     'locus',
     'samplename'], axis=1)
    
    X_array = np.array(strat_train_set.loc[:, strat_train_set.columns != focal_taxon])
    Y_array = np.array(strat_train_set[focal_taxon]).reshape(-1,1)
    
    Xtest_array = np.array(strat_test_set.loc[:, strat_test_set.columns != focal_taxon])
    Ytest_array = np.array(strat_test_set[focal_taxon]).reshape(-1,1)
    
    strat_test_set.shape  
    strat_train_set.shape
    
    ###############################
    # hyperparameter optimization #
    ###############################
    
    #Two options are included here. The first uses
    # Bayesian Sequential Model-based Optimization (SMBO) using HyperOpt. 
    # 
    
    #The second is a random search with scikit
    
    #Note that I am using an easy to use AutoML wrapper for hyperopt called hpsklearn
    #using hyperopt itself is a bit trickier, but not too bad. This way is much
    #simpler though.
    # refit the model
    
    model_random = HyperoptEstimator(regressor=random_forest_regression('random_forest_regression'),
                          algo=tpe.suggest, 
                          max_evals=int(eval_runs), 
                          trial_timeout=30)
    
    # perform the search
    model_random.fit(X=X_array, y=Y_array.ravel())
    
    # FOR RANDOM SEARCH OF HYPERPARAMETERS
    #Code modified from: https://towardsdatascience.com/hyperparameter-tuning-the-random-forest-in-python-using-scikit-learn-28d2aa77dd74
    
    # from sklearn.model_selection import RandomizedSearchCV
    
    # n_estimators = [int(x) for x in np.linspace(start = 200, stop = 500, num = 10)]
    
    # # Number of features to consider at every split
    # max_features = ['auto', 'sqrt']
    
    # # Maximum number of levels in tree
    # max_depth = [int(x) for x in np.linspace(4, 50, num = 2)]
    # max_depth.append(None)
    
    # # Minimum number of samples required to split a node
    # min_samples_split = [2, 5, 10]
    # # Minimum number of samples required at each leaf node
    # min_samples_leaf = [1, 2, 4]
    # # Method of selecting samples for training each tree
    # bootstrap = [True, False]
    # # Create the random grid
    # random_grid = {'n_estimators': n_estimators,
    #                 'max_features': max_features,
    #                 'max_depth': max_depth,
    #                 'min_samples_split': min_samples_split,
    #                 'min_samples_leaf': min_samples_leaf,
    #                 'bootstrap': bootstrap}
    
    # from pprint import pprint
    # pprint(random_grid)
    
    # {'bootstrap': [True, False],
    #                 'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
    #   'max_features': ['auto', 'sqrt'],
    #   'min_samples_leaf': [1, 2, 4],
    #   'min_samples_split': [2, 5, 10],
    #   'n_estimators': [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]}
    
    
    # # define the model
    # model = RandomForestRegressor(
    #                               criterion = 'mse',
    #                               random_state = 666,
    #                               n_jobs = -1
    #                               )
     
    # # FOR RANDOM SEARCH OF HYPERPARAMETERS
    # # # Random search of hyper parameters, using 3 fold cross validation, 
    # # # search across 100 different combinations, and use all available cores
    # model_random = RandomizedSearchCV(estimator = model, 
    #                                 param_distributions = random_grid, 
    #                                 n_iter = 5, #EDIT
    #                                 cv = 3, 
    #                                 verbose=2, 
    #                                 random_state=666, 
    #                                 n_jobs = -1) #-1 means use all processers
    
    #Convert the RandomizedSearchCV object to a randomForestRegressor object
    # tuned_model = RandomForestRegressor(**model_random.best_model())
    
    # #sanity check
    # # yhat = model_random.predict(Xtest_array[:,10:])
    # tuned_model.fit(X=X_array[:,10:], y=Y_array.ravel())
    # yhat2 = tuned_model.predict(Xtest_array[:,10:])
    # r2 = r2_score(Ytest_array,yhat2)
    # r2
    
    ##############
    # model eval #
    ##############
    
    mae = model_random.score(X_array, Y_array.ravel())
    print("MAE: %.3f" % mae)
    
    # summarize the best model
    print(model_random.best_model())
    
    yhat = model_random.predict(Xtest_array)
    
    r2 = r2_score(Ytest_array,yhat)
    r2
    
    #save the model
    
    joblib.dump(model_random, './models/rf_model_object_16s' + focal_taxon)
    
    #To load use this syntax
    #model_random = joblib.load("PATH")
    
    ##############
    #SHAP values
    ############
    
    # Create Tree Explainer object that can calculate shap values
    explainer = shap.TreeExplainer(model_random.best_model()['learner'])
    
    shap_values = explainer.shap_values(X = X_array)
    list_of_labels = list(strat_train_set.loc[:, strat_train_set.columns != focal_taxon].columns)
    
    # shap.summary_plot(shap_values, 
    #                   X_array,
    #                   feature_names=list_of_labels)
    
    # #Extract SHAP values and then remove some features and try rerunning model
    shapValues = pd.DataFrame(shap_values, columns=list_of_labels, 
                              index = trainsamples)
    
    #write SHAP values to disk. 
    shapValues.to_csv(path_or_buf=("./models/shap_values_randomforest_16s" + focal_taxon + ".csv"))
    
    ###############################################
    #Feature removal following SHAP value creation
    ###############################################
    
    #Figure out unimportant features
    duds = shapValues.sum(axis = 0)[abs(shapValues.sum(axis = 0)) < 0.2]
    
    #Here I am searching for taxa in this list of useless features, bc
    #I WANT to keep all taxa
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
    
    #remove a few things, just to be safe. 
    del X, X_array, Y_array, Ytest_array, split, strat_test_set, strat_train_set
    
    #Note that the reason to do this is because if one has a bunch of shit features
    #they will get picked during feature bagging while training the model, thus
    #diluting the effect of the better features. 
    
    X = pd.read_csv("./processedData/imputed_scaled_16S_metadata.csv")
    
    #Remove features that we have deemed not helpful
    for i in l:
        X.drop(i, inplace=True, axis=1)
    
    # Remove diversity from features
    #X.columns.values
    X = X.loc[:,X.columns != 'shannonsISD']
    
    if 'EN' in X.columns:
        X['compartment'] = X['EN'].map({1: 'EN', 0: 'EP'})
        X['sample'] = X['plant.x'].str.cat(X['compartment'], sep='_')
    else:
         X['sample'] = X['plant.x']
         #bc compartment not included have to strip it off the sample field for merge
         taxa['sample'].replace(to_replace="_E[NP]", value="", regex=True, inplace=True)
         
    X_taxa = pd.merge(X, taxa, on='sample')
    
    N = 100
    res = sorted(range(len(X_taxa[focal_taxon])), key = lambda sub: X_taxa[focal_taxon][sub])[-N:]
    
    X_taxa['focal_taxon_onehot'] = 0
    X_taxa.loc[res, 'focal_taxon_onehot'] = 1
    
    #Code from Geron. 
    if 'EN' in X.columns:
        X_taxa['brutal_hack'] = X_taxa["EN"] + X_taxa['focal_taxon_onehot']
        split = StratifiedShuffleSplit(n_splits=1, test_size=(0.3), random_state=(666))
        for train_index, test_index in split.split(X_taxa, X_taxa['brutal_hack']):
            strat_train_set = X_taxa.loc[train_index]
            strat_test_set = X_taxa.loc[test_index]
    else: 
        split = StratifiedShuffleSplit(n_splits=1, test_size=(0.3), random_state=(666))
        for train_index, test_index in split.split(X_taxa, X_taxa['focal_taxon_onehot']):
            strat_train_set = X_taxa.loc[train_index]
            strat_test_set = X_taxa.loc[test_index]
    
    #define response/predictor split
    #First drop various book keeping columns
    testsamples = strat_test_set['samplename']
    trainsamples = strat_train_set['samplename']
    
    strat_test_set = strat_test_set.drop(['focal_taxon_onehot', 'sample',
     'region_site',
     'Unnamed: 0',
     'X',
     'taxon.x',
     'region_site_plant',
     'plant.x',
     'forward_barcode',
     'reverse_barcode',
     'locus',
     'samplename'], axis=1)
    
    strat_train_set = strat_train_set.drop(['focal_taxon_onehot', 'sample',
     'region_site',
     'Unnamed: 0',
     'X',
     'taxon.x',
     'region_site_plant',
     'plant.x',
     'forward_barcode',
     'reverse_barcode',
     'locus',
     'samplename'], axis=1)
    
    if 'brutal_hack' in strat_train_set.columns:
        strat_train_set = strat_train_set.drop(['brutal_hack', 'compartment'], axis=1)
        strat_test_set = strat_test_set.drop(['brutal_hack', 'compartment'], axis=1)
    
    X_array = np.array(strat_train_set.loc[:, strat_train_set.columns != focal_taxon])
    Y_array = np.array(strat_train_set[focal_taxon]).reshape(-1,1)
    
    Xtest_array = np.array(strat_test_set.loc[:, strat_test_set.columns != focal_taxon])
    Ytest_array = np.array(strat_test_set[focal_taxon]).reshape(-1,1)
    
    strat_test_set.shape  
    strat_train_set.shape
    
    # refit the model
    model_random_reducedFeatureSet = HyperoptEstimator(regressor=random_forest_regression('random_forest_regression'),
                              algo=tpe.suggest, 
                              max_evals=int(eval_runs), 
                              trial_timeout=30)
    
    # perform the search
    model_random_reducedFeatureSet.fit(X=X_array, y=Y_array.ravel())
    
    joblib.dump(model_random_reducedFeatureSet, './models/rf_model_object_reducedFeatureSet16s' + focal_taxon)
    
    ##############
    #SHAP values
    ############
    
    # Create Tree Explainer object that can calculate shap values
    explainer = shap.TreeExplainer(model_random_reducedFeatureSet.best_model()['learner'])
    
    shap_values = explainer.shap_values(X = X_array)
    list_of_labels = list(strat_train_set.loc[:, strat_train_set.columns != focal_taxon].columns)
    
    # shap.summary_plot(shap_values, 
    #                   X_array,
    #                   feature_names=list_of_labels)
    
    # #Extract SHAP values and then remove some features and try rerunning model
    shapValues = pd.DataFrame(shap_values, columns=list_of_labels, 
                              index = trainsamples)
    
    #write SHAP values to disk. 
    shapValues.to_csv(path_or_buf=("./models/shap_values_randomforest_reducedFeatureSet16s" + focal_taxon + ".csv"))
    
    #Write the taxon and new r2 to the results file
    yhat = model_random_reducedFeatureSet.predict(Xtest_array)
    
    r2_reduced = r2_score(Ytest_array,yhat)
    r2_reduced
    
    List=[focal_taxon,'full model',r2,'reduced model', r2_reduced]
      
    # Open our existing CSV file in append mode
    # Create a file object for this file
    with open('./models/randomForest_results_16s.csv', 'a') as f_object:
      
        # Pass this file object to csv.writer()
        # and get a writer object
        writer_object = writer(f_object)
      
        # Pass the list as an argument into
        # the writerow()
        writer_object.writerow(List)
      
        #Close the file object
        f_object.close()

