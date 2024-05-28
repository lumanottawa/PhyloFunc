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
# read each distance file
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)
        processed_df = process_csv(file_path)
        result_df = pd.concat([result_df, processed_df],axis=1, ignore_index=False)
#obtain the evaluation results for KNN 
result_df.to_csv('KNN_evaluation_results_normalized.csv', index=False)


