####################################################################
# Author: Shan Qu, Delft University of Technology                ###
# First created: May 2019                                        ###
# This code is propriatary under the Delphi Research Consortium  ###
#                                                                ###
# product: ML workflow for event detection                       ###
# prediction with trained model
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
mat_contents = scipy.io.loadmat('Data/data2_n005.mat')
full_data = mat_contents['data']
mat_contents = scipy.io.loadmat('Data/label2_n005.mat')
label_full_data = mat_contents['label']
mat_contents = scipy.io.loadmat('Data/testdata_1dfeature_n005.mat')
feature_full_data = mat_contents['feature']
mat_contents = scipy.io.loadmat('Data/testdata_2dfeature_n005.mat')
feature_glcm = mat_contents['feature_glcm']

# feature_full_data = np.concatenate((feature_full_data, feature_glcm), 1)
# print feature_full_data.shape

# load operators
with open('Data/clf_no2d.pkl', 'rb') as f:
    clf = pickle.load(f)
with open('Data/selector_no2d.pkl', 'rb') as f:
    selector = pickle.load(f)
with open('Data/rfecv_no2d.pkl', 'rb') as f:
    rfecv = pickle.load(f)

# shuffle the datasets
# feature, label = shuffle(feature_full_data, label_full_data, random_state=0)
feature = feature_full_data
label = label_full_data

# normalizition of features
reg_scaler = preprocessing.StandardScaler().fit(feature)
normalized_feature = reg_scaler.transform(feature)

######################### normailized features / use the whole data to train ###################
feature_test = normalized_feature
label_test = label

######################### feature selection ######################
# univariate feature selection with F test for feature scoring
fselected_feature_test = selector.transform(feature_test)

# rfecv + random forest
fselected_feature_test = rfecv.transform(fselected_feature_test)

######################### start classification #####################

label_pred = clf.predict(fselected_feature_test)
print(classification_report(label_test, label_pred))

scipy.io.savemat('Data/label_pred2_n005_no2d.mat', {'label_pred':label_pred})

# plt.show()

