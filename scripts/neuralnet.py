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
from keras.models import Sequential
from keras.layers import Dense
# from keras.wrappers.scikit_learn import KerasRegressor
# from sklearn.model_selection import cross_val_score
# from sklearn.model_selection import KFold
# from sklearn.preprocessing import StandardScaler
# from sklearn.pipeline import Pipeline

#Load data
data = pd.read_csv('../processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv')

#Explore the data
data.head
data.info()
data.describe() #This is like R's summary command when called on numerical data
data['taxon.x'].value_counts()
#data.columns

data['sla'] = data['area_cm2'] / data['mass_extracted_g']

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
data = data[cols]


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

model.add(Dense(12,input_dim=(1),kernel_initializer='normal',activation='relu'))
#model.add(Dense(8,activation='relu'))
model.add(Dense(1, kernel_initializer='normal'))
model.compile(loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])

#Keras needs Numpy arrays or TensorFlow dataset objects as input. 
#The latter is optimized for really large datasets. Need to learn more about it
#as it could be useful for genomics data. 
X_array = np.array(X.loc[:, X.columns != 'shannonsISD'])
Y_array = np.array(X['shannonsISD']).reshape(1,-1)


#An epoch is one pass through all the rows in the data set
#Batch is how many samples (rows) to consider before updating weights
model.fit(X.iloc[1:2419,], Y_array, epochs=20, batch_size=100)

model.summary()

y_pred= model.predict(np.array(X["area_cm2"]).reshape(1,-1))


import matplotlib.pyplot as plt
X.plot(kind='scatter',
       x='area_cm2',
       y='shannonsISD', title='area vs shannons ISD')
plt.plot(X["area_cm2"], y_pred, color='red', linewidth=3)




#Get working on training dataset
model = Sequential()
#Input dim should be the number of predictors, subtracting one because the response
#is included
#The first integer (12) is the number of nodes in the first dense layer

model.add(Dense(12,input_dim=(1,),kernel_initializer='normal',activation='relu'))
#model.add(Dense(8,activation='relu'))
model.add(Dense(1, kernel_initializer='normal'))
model.compile(loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])

#Keras needs Numpy arrays or TensorFlow dataset objects as input. 
#The latter is optimized for really large datasets. Need to learn more about it
#as it could be useful for genomics data. 
X_array = np.array(X.loc[:, X.columns != 'shannonsISD'])
Y_array = np.array(X['shannonsISD']).reshape(1,-1)


#An epoch is one pass through all the rows in the data set
#Batch is how many samples (rows) to consider before updating weights
model.fit(np.array(X["area_cm2"]).reshape(1,-1), Y_array, epochs=20, batch_size=100)

model.summary()






