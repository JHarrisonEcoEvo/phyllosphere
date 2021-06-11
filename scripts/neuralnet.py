#Code adapted from https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/

import pandas as pd 
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

data = pd.read_csv('processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv')

data.head()


# define the keras model
model = Sequential()

model.add(Dense(12, input_dim=8, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(1, activation='sigmoid'))