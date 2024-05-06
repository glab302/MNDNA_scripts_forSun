import pandas as pd
import numpy as np
import json, pickle, argparse
import math, os, random
from collections import Counter
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, KFold, RepeatedStratifiedKFold, cross_val_score, cross_val_predict, cross_validate
from sklearn.metrics import roc_auc_score
random.seed(1234)
np.random.seed(1234)

filepath = 'Apc_model/data/'
trnval = pd.read_csv(filepath+'2.ApcDataMatrix_SignedFindMarkers_40regions.trnvalset.log',sep='\t'); trnval.index = [i.split('.')[0] for i in trnval.index.tolist()]
test1 = pd.read_csv(filepath+'2.ApcDataMatrix_SignedFindMarkers_40regions.testset.log',sep='\t')
trnval = trnval.T; test1 = test1.T

#### training set
trnval_x = trnval[[i for i in trnval.columns.tolist() if 'A_Apc' in i or 'A_Ctrl' in i]].T
trnval_y = [1 if 'A_Apc' in i else 0 for i in trnval_x.index.tolist()]; trnval_y = pd.DataFrame(trnval_y); trnval_y.index = trnval_x.index.tolist(); trnval_y= trnval_y[0]
savelabel = 'apc13_wt18_'

#### test set
## 0. test
test_rep = test1[[i for i in test1.columns.tolist() if 'A_Apc20' in i or 'A_Ctrl20' in i]].T
## 1. young apc
test_young = pd.concat([trnval[[i for i in trnval.columns.tolist() if 'A_y' in i]].T, test1[[i for i in test1.columns.tolist() if 'A_Apc' in i or 'A_Ctrl' in i]].T])
test_young = test_young.loc[set(test_young.index.tolist())-set(test_rep.index.tolist()),:]
## 2. tranplanted
test_transplanted = pd.concat([trnval[[i for i in trnval.columns.tolist() if 'A_tApc' in i or 'A_tCtrl' in i or '_ot' in i]].T, test1[[i for i in test1.columns.tolist() if 'A_tApc' in i or 'A_tCtrl' in i]].T])
## 3. fetal liver
test_fl = pd.concat([trnval[[i for i in trnval.columns.tolist() if 'A_fl' in i]].T, test1[[i for i in test1.columns.tolist() if 'A_flApc_' in i or 'A_flCtrl_' in i]].T])
## 4. FL+IL18
test_fl18 = test1[[i for i in test1.columns.tolist() if 'A_flApcIl' in i or 'A_flCtrlIl' in i]].T
## 5.1. +PBS
test_PBS = test1[[i for i in test1.columns.tolist() if 'A_PBS' in i or 'W_PBS' in i]].T
## 5.2. +IgG
test_IgG = test1[[i for i in test1.columns.tolist() if 'A_IgG' in i or 'W_IgG' in i]].T
## 5.3. +antiIL18
test_antiIL18 = test1[[i for i in test1.columns.tolist() if 'A_antiIL18' in i or 'W_antiIL18' in i]].T


# model
svc = SVC(kernel='linear', probability=True, random_state=1234)#gamma='auto', 
df_svc = cross_val_score(svc, trnval_x, trnval_y, cv=10, scoring='accuracy').mean()

n_test_set = 0
df_merge = pd.DataFrame()
namelist = ['20w','young','transplanted','20w+PBS','20w+IgG','20w+antiIL18','FetalLiver','FetalLiver+IL18']#
for test_set in [test_rep,test_young,test_transplanted,test_PBS,test_IgG,test_antiIL18,test_fl,test_fl18]:#
    n_splits = 10
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=1234)
    df_train_proba_temp, df_train_pred_temp, df_train_label_temp = [], [], []
    df_valid_proba_temp, df_valid_pred_temp, df_valid_label_temp, df_valid_name_temp = [], [], [], []
    df_test_proba_temp, df_test_pred_temp = pd.DataFrame(), pd.DataFrame()
    # df_test_pred_temp1, df_test_proba_temp1 = pd.DataFrame(), pd.DataFrame()
    round_i = 0
    classifier = svc
    for i, (train, valid) in enumerate(cv.split(trnval_x)):
        print(valid)
        train_X, valid_X, train_y, valid_y = trnval_x.iloc[train],trnval_x.iloc[valid],trnval_y.iloc[train],trnval_y.iloc[valid]
        classifier.fit(train_X,train_y)
        df_train_label_temp.extend(np.array(train_y))
        df_train_pred_temp.extend(classifier.predict(train_X))
        df_train_proba_temp.extend(classifier.predict_proba(train_X)[:,1])
        df_valid_name_temp.extend(valid_X.index.tolist())
        df_valid_label_temp.extend(np.array(valid_y))
        df_valid_pred_temp.extend(classifier.predict(valid_X))
        df_valid_proba_temp.extend(classifier.predict_proba(valid_X)[:,1])
        df_test_pred_temp['test_pred_%i' % (round_i)] = classifier.predict(test_set)
        df_test_proba_temp['test_pred_%i' % (round_i)] = classifier.predict_proba(test_set)[:,1]
        round_i = round_i+1
    df_valid_proba = pd.DataFrame(df_valid_proba_temp); df_valid_proba.index = df_valid_name_temp
    df_valid_proba['testset_label'] = 'cross-validation'
    df_test_proba_temp.index = test_set.index.tolist()
    df_test_proba_mean = pd.DataFrame(df_test_proba_temp.mean(axis=1))
    df_test_proba_mean['testset_label'] = namelist[n_test_set]#'test_pred_%i' % (n_test_set)
    df_merge = pd.concat([df_merge, df_test_proba_mean])
    n_test_set = n_test_set+1


df_merge = pd.concat([df_merge, df_valid_proba])
df_merge['Label'] = 1
ctrl = [i for i in df_merge.index.tolist() if 'Ctrl' in i]
w = [i for i in df_merge.index.tolist() if 'W_' in i]
# w_treat = [i for i in df_merge.index.tolist() if 'A_antiIL18' in i]
df_merge.loc[ctrl,'Label'] = 0
df_merge.loc[w,'Label'] = 0
# df_merge.loc[w_treat,'Label'] = 0
df_merge.to_csv('Apc_model/result/svc_result.log', sep='\t')

roc_auc_score(df_train_label_temp, df_train_proba_temp)
roc_auc_score(df_valid_label_temp, df_valid_proba_temp)

# pd.concat([pd.DataFrame(df_train_label_temp), pd.DataFrame(df_train_proba_temp)], axis=1).to_csv('Apc_model/result/svc_trn.log', sep='\t')
# pd.concat([pd.DataFrame(df_valid_label_temp), pd.DataFrame(df_valid_proba_temp)], axis=1).to_csv('Apc_model/result/svc_val.log', sep='\t')

pd.concat([pd.DataFrame(trnval_x.columns.tolist()), pd.DataFrame(classifier.coef_.tolist()[0])], axis=1).to_csv('Apc_model/result/feature_importance.log',sep='\t')