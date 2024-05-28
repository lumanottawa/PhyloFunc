# Evaluate the performents of four distance matrix by machine learning methods KNN, MLP,SVM
#original data files in the folder 2_data_processing\1_PhyloFunc calculation\3_human gut microbiome\distances:
# 1. KNN classification method for each distance matrix
import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
import os
from sklearn.preprocessing import MinMaxScaler
current_path = os.getcwd()
#get the current path
folder_path = current_path
length=len(folder_path)
#loop the classification 
def process_csv(file_path):
    data = pd.read_csv(file_path, index_col=0)
    X = np.array(data.iloc[:, 0:])
    scaler = MinMaxScaler()
    normalized_data = scaler.fit_transform(X)
    y = np.array([i.split('.')[0] for i in data.index])
    accuracy_drug_list = []
    # use leave-one-out evaluation to assess the accuracy of classification
    loo = LeaveOneOut()
    for k in range(1, 6):
        y_pred_drug_all = []
        y_true_drug_all = []
        knn = KNeighborsClassifier(n_neighbors=k,metric='precomputed',algorithm="brute")
        for train_index, test_index in loo.split(normalized_data):
            X_train, X_test = np.delete(X[train_index],test_index, axis=1),np.delete(X[test_index],test_index,axis=1)
            y_train, y_test = y[train_index], y[test_index]
            
            knn.fit(X_train, y_train)
            y_pred = knn.predict(X_test)
            y_pred_drug = y_pred[0]
            y_true_drug = y_test[0]
            y_pred_drug_all.append(y_pred_drug)
            y_true_drug_all.append(y_true_drug)
        accuracy_drug = accuracy_score(y_true_drug_all, y_pred_drug_all)
        accuracy_drug_list.append(accuracy_drug)
    df = pd.DataFrame({'k': list(range(1,6)), 'accuracy_drug': accuracy_drug_list})
    f=file_path[length+1:]
    df.rename(columns={'accuracy_drug': f}, inplace=True)
    return df

result_df = pd.DataFrame()
# read each distance file in the current path
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)
        processed_df = process_csv(file_path)
        result_df = pd.concat([result_df, processed_df],axis=1, ignore_index=False)
#obtain the evaluation results for KNN 
result_df.to_csv('KNN_evaluation_results_normalized.csv', index=False)

#2. MLP classification method for each distance matrix
import pandas as pd
from sklearn.neural_network import MLPClassifier
import os
from sklearn.model_selection import GridSearchCV, LeaveOneOut
current_path = os.getcwd()
#use name of distance file: Jaccard distance.csv，Euclidean distance.csv, Bray_Curtis.csv to instead PhyloFunc_distance to get the evaluation results
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

#3. SVM classification method for each distance matrix
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.svm import SVC
import os
current_path = os.getcwd()
##use name of distance file: Jaccard distance.csv，Euclidean distance.csv, Bray_Curtis.csv to instead PhyloFunc_distance to get the evaluation results
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
