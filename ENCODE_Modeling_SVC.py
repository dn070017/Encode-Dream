# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 17:56:33 2016

@author: Chester
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import sys
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import SelectKBest, SelectPercentile
from sklearn.feature_selection import chi2
from sklearn.decomposition import PCA
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.cross_validation import KFold
from sklearn import grid_search
from sklearn import linear_model
from sklearn import svm
import sklearn.metrics as skMetric

""""""""""""""""""""""""""""""
# set parameters
""""""""""""""""""""""""""""""
str_filename_inputFile_train = str(sys.argv[1])
str_filename_inputFile_test = str(sys.argv[2])
#str_filename_inputFile_train = "D:/Phd/Grade_04/ENCODE/Data/H1-hESC_ATF2.table_head"
#str_filename_inputFile_train = "C:/Users/Sony/Desktop/ENCODE/Data/Feature/merge_feature_train_ATF2_H1-hESC_head_10000.csv"
#str_filename_inputFile_test = "C:/Users/Sony/Desktop/ENCODE/Data/Feature/merge_feature_test_ATF2_K562_head_10000.csv"

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
warnings.filterwarnings("ignore")
# loading
""""""""""""""""""""""""""""""
#np_header = ['DNase','HelT_X','HelT_Y','HelT_Z','MGW_X','MGW_Y','MGW_Z','ProT_X','ProT_Y','ProT_Z','Roll_X','Roll_Y','Roll_Z','Kmer','GC','Expression']
#np_feature_name = np_header
list_inputFile = []
with open(str_filename_inputFile_train, 'r') as file_inputFile:
    for line in file_inputFile:
        list_line = str(line).strip().split(",")
        list_inputFile.append(list_line)
np_inputFile = np.array(list_inputFile, dtype=np.dtype((str, 10)))
np_feature_train = np_inputFile[:, 3:(np_inputFile.shape[1]-1)]
np_label_train = np_inputFile[:, (np_inputFile.shape[1]-1)]
del np_inputFile
list_inputFile = []
with open(str_filename_inputFile_test, 'r') as file_inputFile:
    for line in file_inputFile:
        list_line = str(line).strip().split(",")
        list_inputFile.append(list_line)
np_inputFile = np.array(list_inputFile, dtype=np.dtype((str, 10)))
np_instanceKey_test = np_inputFile[:, :3]
np_feature_test = np_inputFile[:, 3:np_inputFile.shape[1]]
del np_inputFile
print("01-- loading completed")

# preprocessing
""""""""""""""""""""""""""""""
## count missing value
int_cntMissing = 0
for idx_line in range(0, np_feature_train.shape[0]):
    if "na" in np_feature_train[idx_line,:]:
        int_cntMissing = int_cntMissing + 1
np_feature_train = np_feature_train[np_feature_train[:,1] != 'na',:]
np_label_train = np_label_train[np_feature_train[:,1] != 'na']
print("Missing in training data: " + str(int_cntMissing))

## one hot encoding
#sk_oneHotEncoder = OneHotEncoder()
#np_oneHotEncode = sk_oneHotEncoder.fit_transform(np_feature[:,24].reshape(np_feature[:,24].shape[0],1)).toarray()
#np_feature_name = [list(item) for item in zip(np.repeat("location",np_oneHotEncode.shape[1]), ["{:02d}".format(x) for x in range(np_oneHotEncode.shape[1])])]
#np_feature_name = np.array(["_".join(item) for item in np_feature_name])
#np_feature_name = np.concatenate((np_header[:23], np_header[25:], np_feature_name), axis=0)
#np_feature = np.concatenate((np_feature[:,:23], np_feature[:,25:], np_oneHotEncode), axis=1)
np_feature_train[np_feature_train == "na"] = 0
np_feature_train = np_feature_train.astype(np.dtype('f4'))
np_label_train[np_label_train == "U"] = 0
np_label_train[np_label_train == "A"] = 1
np_label_train[np_label_train == "B"] = 1
np_label_train = np_label_train.astype(float)
np_feature_test[np_feature_test == "na"] = 0
np_feature_test = np_feature_test.astype(np.dtype('f4'))

## missing value imputing (most_frequency)
#sk_missingImputer = Imputer(missing_values=np.nan, strategy='most_frequent')
#np_feature = sk_missingImputer.fit_transform(np_feature)

## min-max scaling
sk_minMaxScaler = MinMaxScaler()
np_feature_train = sk_minMaxScaler.fit_transform(np_feature_train)
np_feature_test = sk_minMaxScaler.transform(np_feature_test)
print("02-- preprocessing completed")

# sample reduction
""""""""""""""""""""""""""""""
np_feature_train = np_feature_train[np_feature_train[:,0] != 0,:]
np_label_train = np_label_train[np_feature_train[:,0]!=0]
print("03-- sample reduction completed")
print("Bind: " + str(np.count_nonzero(np_label_train==0)) + "; Unbind: " + str(len(np_label_train) - np.count_nonzero(np_label_train==0)))

# feature reduction
""""""""""""""""""""""""""""""
## pca analysis before feature selection
#np_color = np_label.astype('|S10')
#for idx in range(0, len(np_color)):
#    if np_color[idx] == '1.0': np_color[idx] = 'b'
#    else: np_color[idx] = 'w'
#pca = PCA()
#reduced_data = pca.fit_transform(np_feature)
#fig = plt.figure(figsize=(10, 10), dpi=100)
#fig.suptitle("PCA before Feature Selection", fontsize=24, fontweight='bold')
#ax = fig.add_subplot(111)
#ax.set_xlabel("Component 1: " + "{:0.2f}".format(pca.explained_variance_ratio_[0]))
#ax.set_ylabel("Component 2: " + "{:0.2f}".format(pca.explained_variance_ratio_[1]))
#ax.scatter(reduced_data[:,0], reduced_data[:,1], c=np_color)
#fig.savefig('PCA_beforeSelection.png')

## feature reduction (univariate using chi-square test)
#np_feature = SelectPercentile(chi2, percentile=100).fit_transform(np_feature, np_label)

## pca analysis after feature selection
#reduced_data = pca.fit_transform(np_feature)
#fig = plt.figure(figsize=(10, 10), dpi=100)
#fig.suptitle("PCA after Feature Selection", fontsize=24, fontweight='bold')
#ax = fig.add_subplot(111)
#ax.set_xlabel("Component 1: " + "{:0.2f}".format(pca.explained_variance_ratio_[0]))
#ax.set_ylabel("Component 2: " + "{:0.2f}".format(pca.explained_variance_ratio_[1]))
#ax.scatter(reduced_data[:,0], reduced_data[:,1], c=np_color)
#fig.savefig('PCA_afterSelection.png')
print("04-- feature reduction completed")

# modeling
""""""""""""""""""""""""""""""
## shuffling
np_feature_train_sparse = coo_matrix(np_feature_train)
np_feature_train, np_feature_train_sparse, np_label_train = shuffle(np_feature_train, np_feature_train_sparse, np_label_train, random_state=0)

cost = [2**x for x in range(0, 1)]
gamma = [2**x for x in range(0, 1)]
parameters = [{'kernel':['linear'], 'C':cost, 'class_weight':['balanced'], 'probability':[True], 'cache_size':[1000]}]
estimator = svm.SVC()
estimator_grid = grid_search.GridSearchCV(estimator, parameters, scoring='f1_weighted')
#estimator_grid.set_params(n_jobs=10, cv=2)
estimator_grid.fit(np_feature_train, np_label_train)
np_predict = estimator_grid.best_estimator_.predict(np_feature_train)
#print(estimator_grid.best_estimator_)
print("Accuracy: " + "{:0.4f}".format(skMetric.accuracy_score(np_label_train, np_predict)) + "; Precision: " + "{:0.4f}".format(skMetric.precision_score(np_label_train, np_predict)) + "; Recall: " + "{:0.4f}".format(skMetric.recall_score(np_label_train, np_predict)))
print("05-- parameters tuning completed")

str_filename_outputFile = "/".join(str_filename_inputFile_test.split("/")[:-1])+ "/" + str_filename_inputFile_test.split("/")[-1] + ".predict"
with open(str_filename_outputFile, 'w') as file_outputFile:
    for idx in range(0,np_instanceKey_test.shape[0]):
        if np_feature_test[idx,0] == 0:
            file_outputFile.write(str(np_instanceKey_test[idx, 0]) + "\t" + str(np_instanceKey_test[idx, 1]) + "\t" + str(np_instanceKey_test[idx, 2]) + "\t" + str(estimator_grid.best_estimator_.predict_proba(np_feature_test[idx,:])[0][0]) + "\n")
        else:
            file_outputFile.write(str(np_instanceKey_test[idx, 0]) + "\t" + str(np_instanceKey_test[idx, 1]) + "\t" + str(np_instanceKey_test[idx, 2]) + "\t" + "0" + "\n")

## k-fold cross validation
kf = KFold(np_feature_train.shape[0], n_folds=2)
## modleing with SVC
list_train_true = []
list_train_predict = []
list_test_true = []
list_test_predict = []
for idx_Tr, idx_Te in kf:
    estimator_svc = svm.SVC(kernel='linear', C=estimator_grid.best_estimator_.C, class_weight='balanced', probability=True, cache_size=1000)
    estimator_svc.fit(np_feature_train[idx_Tr], np_label_train[idx_Tr])
    np_predict = estimator_svc.predict(np_feature_train[idx_Tr])
    for idx_true, idx_predict in zip(list(np_label_train[idx_Tr]), np_predict):
        list_train_true.append(float(idx_true))
        list_train_predict.append(idx_predict)
    np_predict = estimator_svc.predict(np_feature_train[idx_Te])
    for idx_true, idx_predict in zip(list(np_label_train[idx_Te]), np_predict):
        list_test_true.append(float(idx_true))
        list_test_predict.append(idx_predict)

float_accuracy = skMetric.accuracy_score(list_train_true, list_train_predict)
float_precision = skMetric.precision_score(list_train_true, list_train_predict)
float_recall = skMetric.recall_score(list_train_true, list_train_predict)
print("Training Score - Accuracy: " + "{:0.4f}".format(float_accuracy) + "; Precision: " + "{:0.4f}".format(float_precision) + "; Recall: " + "{:0.4f}".format(float_recall))
float_accuracy = skMetric.accuracy_score(list_test_true, list_test_predict)
float_precision = skMetric.precision_score(list_test_true, list_test_predict)
float_recall = skMetric.recall_score(list_test_true, list_test_predict)
print("Testing Score - Accuracy: " + "{:0.4f}".format(float_accuracy) + "; Precision: " + "{:0.4f}".format(float_precision) + "; Recall: " + "{:0.4f}".format(float_recall))
