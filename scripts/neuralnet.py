#Neural net analysis of phyllosphere data
#Code adapted from https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/
#Contents of code:
#Loading and examining data 
#Feature selection
#Feature engineering: dummy variable creation, NA imputation, scaling/centering
#Model definition and compilation
#Checking model performance (standard test/train split)

import pandas as pd
import numpy as np
from tensorflow import keras 
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
# from keras.wrappers.scikit_learn import KerasRegressor
# from sklearn.model_selection import cross_val_score
# from sklearn.model_selection import KFold
# from sklearn.preprocessing import StandardScaler
# from sklearn.pipeline import Pipeline

#Load data
data1 = pd.read_csv('./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv')

#Explore the data
data1.head
data1.info()
data1.describe() #This is like R's summary command when called on numerical data
data1['taxon.x'].value_counts()
#data.columns

data1['sla'] = data1['area_cm2'] / data1['mass_extracted_g']

#Handy code for making histograms
# import matplotlib.pyplot as plt
# data.hist(bins=50,figsize=(20,15))
# plt.show()

#Select the columns that I want to use as features
#Avoiding things that are likely duplicative or boring
#TIP: Note the use of slashes here to improve readability, and allow me to 
#quickly remove certain features.

#I am removing those categorical features that need dummy coding
cols=["area_cm2",\
    "absorbance_420",\
    "absorbance_940",\
    "Ambient_Humidity",\
    "Ambient_Temperature",\
    "B",\
    "circumStem",\
    "compartment.x",\
    "contactless_temp",\
    "deadDown",\
    "densitometer.y",\
    "ecs_initial",\
    "ecs_max",
    "elev_m",\
    "FmPrime",\
    "FoPrime",\
    "Fs",\
    "FvP.FmP",\
    "G",\
    "gH.",\
    "height_sample",\
    "julianDate",\
    "Latitude",\
    "Leaf_Temp_Differential",\
    "LEF",\
    "leaves_extracted",\
    "Light_Intensity..PAR.",\
    "lifehistory",\
    "Longitude",\
    "mass_extracted_g",\
    "MEM1.y",\
    "MEM2.y",\
    "NPQt_MPF",\
    "phenology",\
    "Phi2",\
    "PhiNO",\
    "PhiNPQ",\
    "plant_vol",\
    "pressure",\
    "qL",\
    "R",\
    "Rel_Chl_intensity",\
    "Relative_Chlorophyll",\
    "RFd",\
    "shannons_flora.y",\
    "shrubRich",\
    "sidePlantSampled",\
    "sla",\
    "slope_perc",\
    "SPAD_420_intensity",\
    "SPAD_420",\
    "thickness",\
    "TimeofDay",\
   "taxon.x",\
    "toughness",\
    "treeRich",\
    "waterRetention"\
    # "NPQt",\
    ,"shannonsISD"
    ]

#Just the hot shit
data = data1[cols]


#####################################################################
# Do one hot encoding, conversion to numeric, and scaling/centering #
#####################################################################
#Note that you ca a lot of this inside the model, which is a bit better
#since the model can be ported more easily, since wrangling is internal.

#First we figure out which features are numeric and which are not
num_features = []
categorical_features = []

for i in data:
    if data[i].dtype.kind in 'iufc':
        num_features.extend([i])
    else:
        categorical_features.extend([i])
#The 'iufc' thing means: i int (signed), u unsigned int, f float, c complex. 
  
#Note I explored the pipeline options with sklearn and decided against
#for this use case.
from sklearn.preprocessing import StandardScaler
from sklearn.experimental import enable_iterative_imputer 
from sklearn.impute import IterativeImputer


#Impute missing data, the default model here is a ridge regression
#Note how I have to convert back to a pandas data frame
#
imp = IterativeImputer(max_iter=50, verbose=0)
imp.fit(data[num_features])
imputed_df = imp.transform(data[num_features])

#Scale data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(imputed_df)
imputed_scaled_df = pd.DataFrame(scaled_data, columns=data[num_features].columns)

#Make one-hot encoded categorical variables
taxa = pd.get_dummies(data['taxon.x'])
habit = pd.get_dummies(data['lifehistory'])
compartment = pd.get_dummies(data['compartment.x'])
phenology = pd.get_dummies(data['phenology'])

#Handy way to concatenate data frames, axis decides if fields or rows are added
#can do pretty smart merging, see help.
X = pd.concat([imputed_scaled_df, taxa, habit, compartment, phenology], axis=1)

#Do a final check for Null/NA
Xbool = X.isnull()

for i in Xbool:
    print(Xbool[i].value_counts())

############
# QC step ##
############

#Simulate data that should be predictive of diversity
#This causes model performance to jump up a lot. R2 is .8 roughly, 
#with these setttings.

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

############
# Analysis #
############

#Define the keras model
model = Sequential()
#Input dim should be the number of predictors, subtracting one because the response
#is included
#The first integer (12) is the number of nodes in the first dense layer

#MODEL 1. This model does poorly with our data (R 2 is negative), but does well
#with the simulated data

# model.add(Dense(12,input_dim=strat_train_set.shape[1]-1,activation='relu'))
# model.add(Dense(8,activation='relu'))
# #Linear activation needed for regression
# model.add(Dense(1, activation='linear'))
# model.compile(loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])

#MODEL 2, increasing the number of neurons considerably
model.add(Dense(50,input_dim=strat_train_set.shape[1]-1,activation='relu'))
model.add(Dense(20,activation='relu'))
model.add(Dense(10,activation='relu'))
model.add(Dense(5,activation='relu'))

#Linear activation needed for regression
model.add(Dense(1, activation='linear'))
model.compile(loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])

#Keras needs Numpy arrays or TensorFlow dataset objects as input. 
#The latter is optimized for really large datasets. Need to learn more about it
#as it could be useful for genomics data. 
X_array = np.array(strat_train_set.loc[:, strat_train_set.columns != 'shannonsISD'])
#Note the use of reshape. This is because the array is of dimension [484,]
#The null dimension causes problems, so this reshape makes it have a dimension of [484,1]
#The -1 asks NUmpy to figure out the missing shape...so I could have put [484,1]
#but that would be less robust.
Y_array = np.array(strat_train_set['shannonsISD']).reshape(-1,1)

#An epoch is one pass through all the rows in the data set
#Batch is how many samples (rows) to consider before updating weights
history =  model.fit(X_array, Y_array, epochs=20, batch_size=50)

#Learn a bit about model structure
model.summary()


##############################################
# Assess feature importance with SHAP values. 
# This is prohibitively time consuming 
###############################################
# import shap

# # print the JS visualization code to the notebook
# shap.initjs()
# explainer = shap.KernelExplainer(model = model.predict, data = X_array, link = "identity")

# X_idx = 0

# shap_value_single = explainer.shap_values(X = X_array[X_idx:X_idx+1,:], nsamples = 100)
# X.iloc[X_idx:X_idx+1,:]

# import ipywidgets as widgets
# # Create the list of all labels for the drop down list
# list_of_labels = list(strat_train_set.loc[:, strat_train_set.columns != 'shannonsISD'].columns)

# # Create a list of tuples so that the index of the label is what is returned
# tuple_of_labels = list(zip(list_of_labels, range(len(list_of_labels))))

# # Create a widget for the labels and then display the widget
# current_label = widgets.Dropdown(options=tuple_of_labels,
#                               value=0,
#                               description='Select Label:'
#                               )

# Display the dropdown list (Note: access index value with 'current_label.value')

# print(f'Current Label Shown: {list_of_labels[current_label.value]}\n')

#Nsamples should stay high, but the auto setting is waaay to slow
#Running this for more than a handful of rows is also way to slow. 

# shap_values = explainer.shap_values(X = X_array[0:5,:], nsamples = 100)

# shap.summary_plot(shap_values = shap_values[current_label.value],
#                   features = X_array[0:5,:],
#                   feature_names=list_of_labels
#                   )

# #####################
#Predict to new data#
#####################
X_array_test = np.array(strat_test_set.loc[:, strat_test_set.columns != 'shannonsISD'])
Y_arraytest = np.array(strat_test_set['shannonsISD']).reshape(-1,1)

y_pred_test = model.predict(X_array_test)

#evaluate
# pearsonr(np.squeeze(y_pred_test), np.squeeze(Y_arraytest))

# mean_absolute_error(y_pred_test, Y_arraytest)
# math.sqrt(mean_squared_error(y_pred_test, Y_arraytest))

# plt.scatter(y_pred_test, Y_arraytest)

# r2_score(Y_arraytest,y_pred_test)

# #piss poor

###############

#write the transformed data to disk for future use
    
#Pandas makes writing a lot easier
towrite = pd.concat([data1.iloc[:,0:9], X], axis = 1)
towrite.to_csv(path_or_buf=(".//imputed_scaled_ITS_metadata.csv"))


#########################
# hyperparameter tuning #
#########################

import tensorflow as tf
from tensorflow import keras
import pandas as pd
from tensorflow.keras import layers
from kerastuner.tuners import RandomSearch
import IPython
from tensorflow.keras.models import Sequential
import kerastuner

from tensorflow.keras.layers import Dense
#from tensorflow.keras.optimizers import Adam

def build_model(hp):
    model = keras.Sequential()
    hp_units = hp.Int('units', min_value=10, max_value=100, step=20)
    
    model.add(keras.layers.Dense(units=hp_units, activation='relu'))
    model.add(Dense(1, activation='linear'))

    model.compile(
        #optimizer="adam",
        optimizer=keras.optimizers.Adam(
             hp.Choice('learning_rate', [1e-2, 1e-3, 1e-4])),
             loss='mean_absolute_error',
             metrics=['mean_absolute_error'])

    return model

tuner_band = kerastuner.Hyperband(build_model, objective = 'val_mean_absolute_error',\
                                  max_epochs =2, factor = 3, directory = '.', project_name = 'tuning')


class ClearTrainingOutput(tf.keras.callbacks.Callback):
    def on_train_end(*args, **kwargs):
        IPython.display.clear_output(wait = True)
    
tuner_band.search(X_array, Y_array, epochs = 2, validation_data = (X_array_test, Y_arraytest), \
                  callbacks = [ClearTrainingOutput()])
    
best_hps = tuner_band.get_best_hyperparameters(num_trials = 1)[0]
best_hps.get('learning_rate')

#If not using Jupyter or similar can't display results written in HTML
from IPython.display import display, HTML
display(HTML(tuner_band.results_summary()).data)
