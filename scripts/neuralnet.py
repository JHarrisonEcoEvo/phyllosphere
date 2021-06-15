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

data['sla'] = data['area_cm2'] / data['mass_g']

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
    #"compartment.x",\
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
   # "lifehistory",\
    "Longitude",\
    "mass_extracted_g",\
    "MEM1.y",\
    "MEM2.y",\
    "NPQt_MPF",\
   # "phenology",\
    "Phi2",\
    "PhiNO",\
    "PhiNPQ",\
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
   #"taxon.x",\
    "toughness",\
    "treeRich",\
    "waterRetention"\
    # "NPQt",\
    ]

# split into input (X) and output (Y) variables
X = data[cols]
Y = data['shannonsISD']

#IMPORTANT NOTE: all features need to be numeric, so some sort of encoding
#is needed. Beware of encoding that could be construed as ordinal. This is why
#I do dummy (one hot) encoding.
taxa = pd.get_dummies(data['taxon.x'])
habit = pd.get_dummies(data['lifehistory'])
compartment = pd.get_dummies(data['compartment.x'])
phenology = pd.get_dummies(data['phenology'])

#Handy way to concatenate data frames, axis decides if fields or rows are added
#can do pretty smart merging, see help.
X = pd.concat([X, taxa, habit, compartment, phenology], axis=1)

Need to scale and center and do NAs

#####################
#Do train/test split#
#####################
#Code from Geron. 
#Set seed to get same split. Doing via instantiation of a new pseudorandom
#number generator, instead of changing global seed. Does not matter in this
#case, but this is best practice so you know what rng used for what part of 
#code and will play well with other functions that change seed

np.random.default_rng(666)

def split_train_test(data, test_ratio):
        shuffled_indices = np.random.permutation(len(data))
        test_set_size = int(len(data) * test_ratio)
        test_indices = shuffled_indices[:test_set_size]
        train_indices = shuffled_indices[test_set_size:]
        return data.iloc[train_indices], data.iloc[test_indices]

train_set, test_set = split_train_test(X, 0.2)
print(len(train_set), "train +", len(test_set), "test")

#Incidentally, to do this via scikit learn use this
#from sklearn.model_selection import train_test_split
#train_set, test_set = train_test_split(data, test_size=0.2, random_state=(666))

############
# Analysis #
############

#Define the keras model
model = Sequential()
model.add(Dense(1,input_dim=len(X['area_cm2']),kernel_initializer='normal',activation='relu'))
model.add(Dense(1, kernel_initializer='normal'))
model.compile(loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])

testX = np.asarray(X["area_cm2"]).astype('float32')
testY = np.asarray(Y).astype('float32')

model.fit(testX, testY, epochs=20)
