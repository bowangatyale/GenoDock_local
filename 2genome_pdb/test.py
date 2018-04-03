import os
import sys
import csv
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
#                                              input parameter process                                            #
###################################################################################################################
#1 read command line argument
#Parameters users input:
#1 pdb_id,2 chain_id,3 amino_acid_position,
#4 mutated_amino_acid,5 wildtype_amino_acid,6 bind_site
#7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
#11 germline_somatic
argument_list=sys.argv
pdb_id_u=argument_list[1]
chain_id_u=argument_list[2]
amino_acid_position_u=argument_list[3]
mutated_amino_acid_u=argument_list[4]
wildtype_amino_acid_u=argument_list[5]
bind_site_u=argument_list[6]
SIFT_score_u=argument_list[7]
PPH_score_u=argument_list[8]
GERP_score_u=argument_list[9]
allele_freq_u=argument_list[10]
germline_somatic_u=argument_list[11]


###################################################################################################################
#                                           parse_user_features                                                   #
###################################################################################################################
#Randomforest model required features:
#1 SIFT_score,2 PPH_score,3 gerp_score,
#4 volumn_change_index,5 polarity_change_index
#6 allele_freq,7 bind_site_F,8 bind_site_T,
#9 germline_somatic_g,10 germline_somatic_s

#User Input
#1 pdb_id,2 chain_id,3 amino_acid_position,4 mutated_amino_acid,
#5 wildtype_amino_acid,6 bind_site,7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
#11 germline_somatic
def parse_user_input(mutated_amino_acid,wildtype_amino_acid,
                     bind_site,SIFT_score,PPH_score,
                     gerp_score,allele_freq,germline_somatic):
    #1 compute user polarity_change_index_u
    aa_count=pd.read_excel('AA_count.xlsx').drop(20)
    acid_change=np.array([mutated_amino_acid_u,wildtype_amino_acid_u]).reshape(1,-1)
    acid_change=pd.DataFrame(acid_change,
                            columns=[['mutated_amino_acid','wildtype_amino_acid']])
    aa_count1=aa_count[['AA_name','Polarity_Index']]
    aa_count1.columns=['mutated_amino_acid','mutated_polarity_index']
    aa_count2=aa_count[['AA_name','Polarity_Index']]
    aa_count2.columns=['wildtype_amino_acid','wildtype_polarity_index']

    polarity_merge_result=acid_change.merge(aa_count1,how='left',on=['mutated_amino_acid']).\
    merge(aa_count2,how='left',on=['wildtype_amino_acid'])

    mutated_polarity_index_u=polarity_merge_result['mutated_polarity_index']
    wildtype_polarity_index_u=polarity_merge_result['wildtype_polarity_index']
    polarity_change_index_u=float(mutated_polarity_index_u-wildtype_polarity_index_u)

    #2 compute user volume_change_index_u
    aa_count1=aa_count[['AA_name','VDW_Volume']]
    aa_count1.columns=['mutated_amino_acid','mutated_volume']

    aa_count2=aa_count[['AA_name','VDW_Volume']]
    aa_count2.columns=['wildtype_amino_acid','wildtype_volume']

    volume_merge_result=acid_change.merge(aa_count1,how='left',on=['mutated_amino_acid']).\
    merge(aa_count2,how='left',on=['wildtype_amino_acid'])

    mutated_volume_u=volume_merge_result['mutated_volume']
    wildtype_volume_u=volume_merge_result['wildtype_volume']
    volume_change_index_u=float(np.log2(float(mutated_volume_u)/wildtype_volume_u))

    #3 encode user bind_site
    if bind_site=='True':
        bind_site_T=1
        bind_site_F=0
    else:
        bind_site_T=0
        bind_site_F=1

    #4 encode user germline_somatic
    if germline_somatic=='germline':
        germline_somatic_g=1
        germline_somatic_s=0
    else:
        germline_somatic_g=0
        germline_somatic_s=1
    #5 construct user feature
    user_features=pd.Series([SIFT_score, PPH_score, gerp_score,\
                             volume_change_index_u,polarity_change_index_u,\
                             allele_freq,bind_site_F, bind_site_T,\
                             germline_somatic_g, germline_somatic_s]).values.reshape(1,-1)
    user_features=pd.DataFrame(user_features)

    training_features_name=['SIFT_score', 'PPH_score', 'GERP_score',\
                             'volume_change_index','polarity_change_index',\
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
    global pdb_id_u
    global chain_id_u
    global amino_acid_position_u
    global mutated_amino_acid_u
    global wildtype_amino_acid_u
    global bind_site_u
    global SIFT_score_u
    global PPH_score_u
    global GERP_score_u
    global allele_freq_u
    global germline_somatic_u


    ##1 Check pdb_id
    if (len(pdb_id_u)!=4):
        print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('pdb_id'))
        sys.exit(1)
    ##2 Check chain_id
    if (len(chain_id_u)!=1 or chain_id_u.isupper()=='False'):
        print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('chain_id'))
        sys.exit(1)
    ##3 Check amino_acid_position
    if amino_acid_position_u.isnumeric():
        amino_acid_position_u=int(amino_acid_position_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('amino_acid_position'))
        sys.exit(1)
    ##4 Check mutated_amino_acid
    amino_acid_list=['ALA','LYS','GLU','ASP','ARG','TYR','TRP','THR','SER','HIS',
                     'GLN','CYS','ASN','VAL','PRO','PHE','MET','LEU','ILE','GLY']
    if (mutated_amino_acid_u not in amino_acid_list):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('mutated_amino_acid'))
        sys.exit(1)
    ##5 Check wildtype_amino_acid
    if (wildtype_amino_acid_u not in amino_acid_list):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('wildtype_amino_acid'))
        sys.exit(1)
        sys.exit(1)
    ##6 Check bind_site
    if (bind_site_u=='True') or (bind_site_u=='False'):
        bind_site_u=bool(bind_site_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('bind_site'))
        sys.exit(1)
    ##7 Check SIFT_score
    if SIFT_score_u.replace('.', '', 1).isnumeric():
        SIFT_score_u=float(SIFT_score_u)
        if SIFT_score_u<0 or SIFT_score_u>1:
            print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('SIFT_score'))
            sys.exit(1)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('SIFT_score'))
        sys.exit(1)
    ##8 Check PPH_score
    if PPH_score_u.replace('.', '', 1).isnumeric():
        PPH_score_u=float(PPH_score_u)
        if PPH_score_u<0 or PPH_score_u>1:
            print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit'.format('PPH_score'))
            sys.exit(1)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('PPH_score'))
        sys.exit(1)
    ##9 Check Gerp_score
    if GERP_score_u.replace('.', '', 1).isnumeric():
        GERP_score_u=float(GERP_score_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('GERP_score'))
        sys.exit(1)

    ##10 Check allele_freq
    if allele_freq_u.replace('.', '', 1).isnumeric():
        allele_freq_u=float(allele_freq_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('allele_freq'))
        sys.exit(1)
    ##11 Check germline_somatic
    if (germline_somatic_u != 'germline' and germline_somatic_u != 'somatic'):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('germline_somatic'))
        sys.exit(1)
    ###################################################################################################################
    #                                           finish  input_parameter_check                                         #
    ###################################################################################################################


    clf= joblib.load("trained_random_forest.m")
    user_features,training_features_name=parse_user_input(mutated_amino_acid_u,wildtype_amino_acid_u,
                                                          bind_site_u,SIFT_score_u,PPH_score_u,
                                                          GERP_score_u,allele_freq_u,germline_somatic_u)

    pred_proba=clf.predict_proba(np.array(user_features).reshape(1,-1))

    print ('Random Forest Features:')
    for name, value in zip(training_features_name,user_features.values.tolist()[0]):
        print (name,':',value)

    print ('')
    print ('Probability that binding_affinity_score_change_class is positive: ',round(float(pred_proba[:,1]),3)*100,'%')
    print ('Probability that binding_affinity_score_change_class is negative or zero: ',round(float(pred_proba[:,0]),3)*100,'%')


main()
