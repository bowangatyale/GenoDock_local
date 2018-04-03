import os
import sys
import pandas as pd
import numpy as np
from IPython.display import display
from scipy import stats, integrate
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

#hide warning
pd.options.mode.chained_assignment = None  # default='warn'
###################################################################################################################
#                                          random forest training                                                 #
###################################################################################################################
def random_forest_training():
    #1read csv
    df=pd.read_csv('all_exac_tcga.csv')

    #2create mutation_change_class feature
    binding_affinity_score_change_class=[]
    for ba in df['BA_Change_Vina'].values:
        if ba >0:
            binding_affinity_score_change_class.append('pos')
        else:
            binding_affinity_score_change_class.append('zero_and_neg')
            df.loc[:,('binding_affinity_score_change_class')]=pd.Series(binding_affinity_score_change_class)
    #3 set model features and target
    features=df[['SIFT_score','PPH_score','gerp_score','bind_site',\
                 'volume_change_index','polarity_change_index',\
                 'allele_freq','germline_somatic']]
    target=df['binding_affinity_score_change_class']
    #4.1 encode the target to 0,1
    target_encoded=[]
    for ba in target:
            if ba=='pos':
                target_encoded.append(1)
            else:
                target_encoded.append(0)
    #4.2 encode bind_site to bind_site_T and bind_site_F
    bind_site_class=[]
    for value in features['bind_site']:
        if value==True:
            bind_site_class.append('T')
        else:
            bind_site_class.append('F')
    features_encoded=features
    features_encoded.loc[:,('bind_site')]=bind_site_class
    #4.3 encode germline_somatic to germline_somatic_g and germline_somatic_s
    germline_somatic_class=[]
    for value in features['germline_somatic']:
        if value=='germline':
            germline_somatic_class.append('g')
        elif value=='somatic':
            germline_somatic_class.append('s')
    features_encoded.loc[:,('germline_somatic')]=germline_somatic_class

    features_encoded=pd.get_dummies(features_encoded)  #encode

    #5 train the random_forest model with all data
    clf=RandomForestClassifier(random_state=10000,n_estimators=1500,max_features='log2',min_samples_leaf=30)
    clf=clf.fit(features_encoded,target_encoded)
    return clf


print ('Start Train Random Forest Model.')
clf=random_forest_training()
print ('Finish Train Random Forest Model.')

joblib.dump(clf, "trained_random_forest.m")
