#MLP classification method for each distance matrix
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neural_network import MLPClassifier
import os
current_path = os.getcwd()
file_path ="".join([current_path,"/PhyloFunc_distance.csv"])
np.random.seed(1)
data = pd.read_csv(file_path, index_col=0)
# normalize distances
X = np.array(data.iloc[:, 0:])
scaler = MinMaxScaler()
normalized_data = scaler.fit_transform(X)
y = np.array([i.split('.')[0] for i in data.index])
# construct models
mlp = MLPClassifier()
hidden_layer_sizes_range = [(layer1, layer2) for layer1 in range(4, 12, 2) for layer2 in range(3, 8, 2)]
# set parameters
param_grid = {
    'hidden_layer_sizes': hidden_layer_sizes_range,
    'activation': ['relu', 'tanh'],
    'solver': ['adam', 'sgd','lbfgs'],
    'alpha': [0.001, 0.01],
    'learning_rate': ['constant', 'adaptive'],
}
loo = LeaveOneOut()
# search paramaters with GridSearchCV
grid_search = GridSearchCV(mlp, param_grid, cv=loo)
grid_search.fit(normalized_data, y)

# output the best parameters and accuracy
print("Best parameters found: ", grid_search.best_params_)
print("Best accuracy: ", grid_search.best_score_)
#use name of distance file: Jaccard distance.csv，Euclidean distance.csv, Bray_Curtis.csv to instead PhyloFunc_distance to get the evaluation results
