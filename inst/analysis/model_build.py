import sys
from tensorflow import keras as K
import tensorflow as tf
from tensorflow.keras import regularizers
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
import matplotlib.pyplot as plt
from tensorflow.keras.backend import clear_session
from sklearn.model_selection import PredefinedSplit
import math
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV

data = pd.read_csv("./model_data213/GSE673_pathway_score.csv")


print(data["type"].value_counts())
data.type = data.type.astype(str).map({'malignant': 0, 'normal': 1})
x = data.drop('type', axis=1)
y = data['type']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=.2, random_state=10)
feature = list(x.columns)


# build model
clear_session()
model = K.models.Sequential()
model.add(K.layers.Dense(units=300, input_dim=213, activation='sigmoid'))
model.add(K.layers.Dropout(0.3))
model.add(K.layers.Dense(units=200, activation='sigmoid'))
model.add(K.layers.Dropout(0.2))
model.add(K.layers.Dense(units=100, activation='sigmoid'))
model.add(K.layers.Dropout(0.1))
model.add(K.layers.Dense(units=10, activation='sigmoid'))
model.add(K.layers.Dense(units=1, activation='sigmoid'))
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

b_size = 50
max_epochs = 100
h = model.fit(x_train, y_train, batch_size=b_size, epochs=max_epochs, shuffle=True, verbose=1)  ##调参，节点数，层数
eval = model.evaluate(x_train, y_train, verbose=0, batch_size=b_size)
print("Evaluation on train data: loss = %0.6f accuracy = %0.2f%% \n" % (eval[0], eval[1] * 100))

eval = model.evaluate(x_test, y_test, verbose=0, batch_size=b_size)
print("Evaluation on test data: loss = %0.6f accuracy = %0.2f%% \n" % (eval[0], eval[1] * 100))

