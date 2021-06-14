#Code adapted from https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/

import pandas as pd
import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

data = pd.read_csv('processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv')

data.head()

#Select the columns that I want to use as features in the net
cols=["taxon.x","area_cm2","mass_extracted_g","circumStem","width","plantMeasurements__width2","plant_vol","plantMeasurements__height","height_sample","Ambient_Humidity","Ambient_Temperature","Leaf_Temp_Differential","LEF","Light_Intensity..PAR.","NPQt","Phi2","PhiNO","PhiNPQ","Relative_Chlorophyll","thickness","absorbance_420","absorbance_530","absorbance_605","absorbance_650","absorbance_730","absorbance_850","absorbance_880","absorbance_940","contactless_temp","ecs_initial","ecs_max","ecs_r_squared","FmPrime","FoPrime","Fs","FvP.FmP","G","gH.","NPQt_MPF","PhiNPQ","pressure","qL","R","Rel_Chl_intensity","RFd","SPAD_420","SPAD_420_intensity","SPAD_530","SPAD_530_intensity","SPAD_605","SPAD_605_intensity","SPAD_650","SPAD_650_intensity","SPAD_730","SPAD_730_intensity","SPAD_850","SPAD_850_intensity","SPAD_880","SPAD_880_intensity","spad_raw1","spad_raw2","spad_raw3","TimeofDay","Latitude","Longitude","waterRetention","toughness","elev_m","slope_perc","treeRich","shrubRich","deadDown","latitude","altitude","longitude","densitometer.y","shannons_flora.y","MEM1.y","MEM2.y","julianDate"]

# split into input (X) and output (Y) variables
X = data[cols]
Y = data['shannonsISD']

#Make dummy variables (one hot encoding) for all non-numeric variables in the data
taxa = pd.get_dummies(X['taxon.x'])
X = [X, taxa]
X = pd.concat(X, axis=1)

#lots of NULLs
#print(X[X.isnull().any(axis=1)])

# define the keras model
	# create model
model = Sequential()
model.add(Dense(1,input_dim=len(X['area_cm2']),kernel_initializer='normal',activation='relu'))
model.add(Dense(1, kernel_initializer='normal'))
model.compile(loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])

testX = np.asarray(X["area_cm2"]).astype('float32')
testY = np.asarray(Y).astype('float32')

model.fit(testX, testY, epochs=20)

