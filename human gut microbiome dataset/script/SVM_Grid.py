#SVM classification method for each distance matrix
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.svm import SVC
import os
current_path = os.getcwd()
#read Euclidean distance, evaluating other distance only change the file_path
file_path ="".join([current_path,"/PhyloFunc_distance.csv"])
data = pd.read_csv(file_path, index_col=0)
#normalize data
X = np.array(data.iloc[:, 0:])
scaler = MinMaxScaler()
normalized_data = scaler.fit_transform(X)
#get the character y
y = np.array([i.split('.')[0] for i in data.index])
# construct the classification model
svm = SVC()
#set the parameters
param_grid = {
    'C': [0.1, 1, 10,100,1000],
    'gamma': ['auto',0.1, 0.01, 0.001,1,2],
    'kernel': ['linear','rbf','poly','sigmoid','precomputed'],
}
#leave one out evaluation
loo = LeaveOneOut()
# use GridSearchCV to search parameters
grid_search = GridSearchCV(svm, param_grid, cv=loo, scoring='accuracy')
grid_search.fit(normalized_data, y)

# output the best parameters and accuracy
print("Best parameters found: ", grid_search.best_params_)
print("Best accuracy: ", grid_search.best_score_)
#use name of distance file: Jaccard distance.csv，Euclidean distance.csv, Bray_Curtis.csv to instead PhyloFunc_distance to get the evaluation results
