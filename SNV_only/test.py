import os
import sys
import pandas as pd
import numpy as np
from IPython.display import display
from scipy import stats, integrate
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
#hide warning
pd.options.mode.chained_assignment = None  # default='warn'

from sklearn.externals import joblib

###################################################################################################################
#                                              input parameter process                                            #
###################################################################################################################
#1 read command line argument
#Parameters users input:
#1 SIFT_score,2 PPH_score,3 GERP_score
#4 bind_site,5 allele_freq,6 germline_somatic

argument_list=sys.argv
SIFT_score_u=argument_list[1]
PPH_score_u=argument_list[2]
GERP_score_u=argument_list[3]
bind_site_u=argument_list[4]
allele_freq_u=argument_list[5]
germline_somatic_u=argument_list[6]

###################################################################################################################
#                                           parse_user_features                                                   #
###################################################################################################################
#Randomforest model required features:
#1 SIFT_score,2 PPH_score,3 gerp_score,
#4 volumn_change_index,5 polarity_change_index
#6 allele_freq,7 bind_site_F,8 bind_site_T,
#9 germline_somatic_g,10 germline_somatic_s

#User Input
#1SIFT_score，2PPH_score,3GERP_score,
#4bind_site,5allele_freq,6germline_somatic
def parse_user_input(SIFT_score,PPH_score,gerp_score,\
                     bind_site,allele_freq,germline_somatic):

    #1 encode user bind_site
    if bind_site=='True':
        bind_site_T=1
        bind_site_F=0
    else:
        bind_site_T=0
        bind_site_F=1

    #2 encode user germline_somatic
    if germline_somatic=='germline':
        germline_somatic_g=1
        germline_somatic_s=0
    else:
        germline_somatic_g=0
        germline_somatic_s=1
    #3 construct user feature
    user_features=pd.Series([SIFT_score, PPH_score, gerp_score,\
                             allele_freq,bind_site_F, bind_site_T,\
                             germline_somatic_g, germline_somatic_s]).values.reshape(1,-1)
    user_features=pd.DataFrame(user_features)

    training_features_name=['SIFT_score', 'PPH_score', 'GERP_score',\
                            'allele_freq','bind_site_F', 'bind_site_T',\
                            'germline_somatic_g', 'germline_somatic_s']


    return user_features,training_features_name
###################################################################################################################
#                                                   main                                                          #
###################################################################################################################
def main():
    ###################################################################################################################
    #                                                   input_parameter_check                                         #
    ###################################################################################################################
    global SIFT_score_u
    global PPH_score_u
    global GERP_score_u
    global bind_site_u
    global allele_freq_u
    global germline_somatic_u
    ##1 Check SIFT_score
    if SIFT_score_u.replace('.', '', 1).isnumeric():
        SIFT_score_u=float(SIFT_score_u)
        if SIFT_score_u<0 or SIFT_score_u>1:
            print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('SIFT_score'))
            sys.exit(1)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('SIFT_score'))
        sys.exit(1)
    ##2 Check PPH_score
    if PPH_score_u.replace('.', '', 1).isnumeric():
        PPH_score_u=float(PPH_score_u)
        if PPH_score_u<0 or PPH_score_u>1:
            print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit'.format('PPH_score'))
            sys.exit(1)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('PPH_score'))
        sys.exit(1)
    ##3 Check Gerp_score
    if GERP_score_u.replace('.', '', 1).isnumeric():
        GERP_score_u=float(GERP_score_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('GERP_score'))
        sys.exit(1)
    ##4 Check bind_site
    if (bind_site_u=='True') or (bind_site_u=='False'):
        bind_site_u=bool(bind_site_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('bind_site'))
        sys.exit(1)
    ##5 Check allele_freq
    if allele_freq_u.replace('.', '', 1).isnumeric():
        allele_freq_u=float(allele_freq_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('allele_freq'))
        sys.exit(1)
    ##6 Check germline_somatic
    if (germline_somatic_u != 'germline' and germline_somatic_u != 'somatic'):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('germline_somatic'))
        sys.exit(1)
    ###################################################################################################################
    #                                           finish  input_parameter_check                                         #
    ###################################################################################################################
    clf = joblib.load("trained_random_forest.m")

    user_features,training_features_name=parse_user_input(SIFT_score_u,PPH_score_u,GERP_score_u,
                                                          bind_site_u,allele_freq_u,germline_somatic_u)

    pred_proba=clf.predict_proba(np.array(user_features).reshape(1,-1))

    print ('Random Forest Features:')
    for name, value in zip(training_features_name,user_features.values.tolist()[0]):
        print (name,':',value)

    print ('')
    print ('Probability that binding_affinity_score_change_class is positive: ',round(float(pred_proba[:,1]),3)*100,'%')
    print ('Probability that binding_affinity_score_change_class is negative or zero: ',round(float(pred_proba[:,0]),3)*100,'%')

main()
