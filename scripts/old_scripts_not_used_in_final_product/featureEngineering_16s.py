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
data1 = pd.read_csv('../processedData/16smetadat_wrangled_for_post_modeling_analysis.csv')

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
    "compartment.y",\
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
   "taxon_final",\
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
#Note that you can do a lot of this inside the model, which is a bit better
#since the model can be ported more easily, since wrangling is internal. I am 
#not doing this for the time being.

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
taxa = pd.get_dummies(data['taxon_final'])
habit = pd.get_dummies(data['lifehistory'])
compartment = pd.get_dummies(data['compartment.y'])
phenology = pd.get_dummies(data['phenology'])

#Handy way to concatenate data frames, axis decides if fields or rows are added
#can do pretty smart merging, see help.
X = pd.concat([imputed_scaled_df, taxa, habit, compartment, phenology], axis=1)

#Do a final check for Null/NA
Xbool = X.isnull()

for i in Xbool:
    print(Xbool[i].value_counts())

#write the transformed data to disk for future use
    
#Pandas makes writing a lot easier
towrite = pd.concat([data1.iloc[:,0:9], X], axis = 1)
towrite.to_csv(path_or_buf=("../processedData/imputed_scaled_16S_metadata.csv"))

