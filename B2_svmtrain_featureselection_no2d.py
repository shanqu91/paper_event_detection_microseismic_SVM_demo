####################################################################
# Author: Shan Qu, Delft University of Technology                ###
# First created: May 2019                                        ###
# This code is propriatary under the Delphi Research Consortium  ###
#                                                                ###
# product: ML workflow for event detection                       ###
# training with only 1D features
####################################################################


import pandas as pd
import numpy as np
from sklearn import svm, datasets, preprocessing, feature_selection, ensemble
import pickle
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.utils import shuffle
import librosa
import scipy
import math
import cmath
from statsmodels import robust

# load full data and label and features
mat_contents = scipy.io.loadmat('Data/data1_n005.mat')
full_data = mat_contents['data']
mat_contents = scipy.io.loadmat('Data/label1_n005.mat')
label_full_data = mat_contents['label']
mat_contents = scipy.io.loadmat('Data/traindata_1dfeature.mat')
feature_full_data = mat_contents['feature']
mat_contents = scipy.io.loadmat('Data/traindata_2dfeature.mat')
feature_glcm = mat_contents['feature_glcm']

# feature_full_data = np.concatenate((feature_full_data, feature_glcm), 1)
# print feature_full_data.shape

# shuffle the datasets
feature, label = shuffle(feature_full_data, label_full_data, random_state=0)

# normalizition of features
normalized_feature = preprocessing.StandardScaler().fit_transform(feature)

######################### normailized features / use the whole data to train ###################
# feature_train, feature_test, label_train, label_test = train_test_split(normalized_feature, label, test_size=0.25, random_state=0)
feature_train = normalized_feature
label_train = label

######################### feature selection ######################
# univariate feature selection with F test for feature scoring
selector = feature_selection.SelectPercentile(feature_selection.f_classif, percentile=30)
selector = selector.fit(feature_train, label_train)
fselected_feature_train = selector.transform(feature_train)

# rfecv + random forest
rfecv = feature_selection.RFECV(ensemble.RandomForestClassifier(), step=1, cv=5, scoring='accuracy')
rfecv = rfecv.fit(fselected_feature_train, label_train)
fselected_feature_train = rfecv.transform(fselected_feature_train)

print("Optimal number of features : %d" % rfecv.n_features_)

# # Plot number of features VS. cross-validation scores
# plt.figure()
# plt.xlabel("Number of features selected")
# plt.ylabel("Cross validation score")
# plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)

# plt.show()

######################### start training using SVM ######################

parameter = [{'kernel': ['rbf'],
              'gamma': ['auto'],
              'C': np.logspace(-3, 3, 10),
              'class_weight': ['balanced']
              }]

clf = GridSearchCV(svm.SVC(decision_function_shape='ovr'), parameter, cv=5)
clf = clf.fit(fselected_feature_train, label_train)

print("The best parameters are %s with a score of %0.2f"
      % (clf.best_params_, clf.best_score_))

means = clf.cv_results_['mean_test_score']
stds = clf.cv_results_['std_test_score']
# print means
# print stds

# label_true, label_pred = label_test, clf.predict(selector.transform(feature_test))
# print(classification_report(label_true, label_pred))


with open('Data/clf_no2d.pkl', 'wb') as f:
    pickle.dump(clf, f)
with open('Data/rfecv_no2d.pkl', 'wb') as f:
    pickle.dump(rfecv, f)
with open('Data/selector_no2d.pkl', 'wb') as f:
    pickle.dump(selector, f)

plt.show()
