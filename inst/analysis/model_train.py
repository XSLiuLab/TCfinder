#########################################################  加载模块
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
import pandas as pd
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
########################################################################### Build partition function


def FindLayerNodesLinear(n_layers, first_layer_nodes, last_layer_nodes):
    layers = []
    nodes_increment = (last_layer_nodes - first_layer_nodes)/ (n_layers-1)
    nodes = first_layer_nodes
    for i in range(1, n_layers+1):
        layers.append(math.ceil(nodes))
        nodes = nodes + nodes_increment
    return layers

def FinddropoutLinear(n_layers, dropout):
    layers = []
    nodes_increment = round(dropout/(n_layers-1),2)
    nodes = dropout
    for i in range(1, n_layers+1):
        layers.append(nodes)
        nodes = round(nodes - nodes_increment,2)
        if(nodes <= 0):
            nodes = 0
    
    return layers


def createmodel(n_layers, first_layer_nodes, last_layer_nodes, activation_func, loss_func,dropout):
    model = Sequential()
    n_nodes = FindLayerNodesLinear(n_layers, first_layer_nodes, last_layer_nodes)
    n_dropout = FinddropoutLinear(n_layers,dropout)
    for i in range(1, n_layers):
        if i==1:
            model.add(Dense(first_layer_nodes, input_dim=train_x.shape[1], activation=activation_func))
            model.add(K.layers.Dropout(rate=n_dropout[0]))
        else:
            model.add(Dense(n_nodes[i-1], activation=activation_func))
            model.add(K.layers.Dropout(rate=n_dropout[i-1]))
    model.add(Dense(train_y.shape[1], activation='softmax'))
    model.compile(optimizer='adam', loss=loss_func, metrics = ["accuracy"]) #note: metrics could also be 'mse'
    
    return model


train_x = pd.read_csv("/public/slst/home/wuchx/project/mcIdentify/mcIdentify/code/train673_model/train_x.csv")
train_y = pd.read_csv("/public/slst/home/wuchx/project/mcIdentify/mcIdentify/code/train673_model/train_x.csv")
train_x1, value_x, train_y1, value_y = train_test_split(train_x, train_y,train_size=0.8, test_size=0.2, random_state=1)


############################################################################ build function
train_val_features = np.concatenate((train_x,value_x),axis=0)
train_val_labels = np.concatenate((train_y,value_y),axis=0)
test_fold = np.zeros(train_val_features.shape[0]) 
test_fold[:train_x1.shape[0]] = -1  
ps = PredefinedSplit(test_fold=test_fold)
#####################################################################Set parameter range
model =  KerasClassifier(build_fn=createmodel, verbose = False)  

activation_funcs = ['sigmoid', 'relu'] 
#activation_funcs = ['relu'] 
loss_funcs = ['binary_crossentropy']

param_grid = dict(n_layers=[3,4,5,6], first_layer_nodes = [200,250,300,350,400,450,500], last_layer_nodes = [10,20,30], dropout=[0.1,0.3,0.5], activation_func = activation_funcs, loss_func = loss_funcs, batch_size = [100,80,50], epochs = [50,100])

grid = GridSearchCV(estimator = model, param_grid = param_grid,cv=ps,n_jobs=1)
#grid = RandomizedSearchCV (estimator = model, param_grid = param_grid,cv=3,n_jobs=5)
################################################################ trainning
grid.fit(train_val_features, train_val_labels)
############################################################### Output the highest accuracy and corresponding parameters
print(grid.best_score_)
print(grid.best_params_)
f = open("/public/slst/home/wuchx/project/mcIdentify/mcIdentify/code/train673_model/model_result.txt") 
f.write(str(grid.best_score_))
f.write("\n")
f.write(str(grid.best_params_)) 
f.write("\n") 
f.close()

