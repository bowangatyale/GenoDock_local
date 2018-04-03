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
#                                              input parameter read                                               #
###################################################################################################################
#Parameters users input:
#1 mutated_amino_acid,2 wildtype_amino_acid
#3 SIFT_score,4 PPH_score,5 GERP_score,6 allele_freq,7 bind_site
#8 germline_somatic,9 molecular_weight,10 H_bond_donor
#11 H_bond_acceptor,12 rotatable_bond,13 Polar_Surface_Area
argument_list=sys.argv
mutated_amino_acid_u=argument_list[1]
wildtype_amino_acid_u=argument_list[2]
SIFT_score_u=argument_list[3]
PPH_score_u=argument_list[4]
GERP_score_u=argument_list[5]
allele_freq_u=argument_list[6]
bind_site_u=argument_list[7]
germline_somatic_u=argument_list[8]

molecular_weight_u=argument_list[9]
H_bond_donor_u=argument_list[10]
H_bond_acceptor_u=argument_list[11]
rotatable_bond_u=argument_list[12]
Polar_Surface_Area_u=argument_list[13]


###################################################################################################################
#                                           parse_user_features                                                   #
###################################################################################################################
#Randomforest model required features:
#1 SIFT_score,2 PPH_score,3 gerp_score,
#4 volumn_change_index,5 polarity_change_index
#6 allele_freq,7 bind_site_F,8 bind_site_T,
#9 germline_somatic_g,10 germline_somatic_s

#User Input
#1 mutated_amino_acid,2 wildtype_amino_acid
#3 SIFT_score,4 PPH_score,5 GERP_score,6 allele_freq,7 bind_site
#8 germline_somatic,9 molecular_weight,10 H_bond_donor
#11 H_bond_acceptor,12 rotatable_bond,13 Polar_Surface_Area
def parse_user_input(mutated_amino_acid_u,wildtype_amino_acid_u,\
                    SIFT_score_u,PPH_score_u,GERP_score_u,allele_freq_u,\
                    bind_site_u,germline_somatic_u,\
                    molecular_weight_u,H_bond_donor_u,H_bond_acceptor_u,rotatable_bond_u,Polar_Surface_Area_u):
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
    if bind_site_u=='True':
        bind_site_u_T=1
        bind_site_u_F=0
    else:
        bind_site_u_T=0
        bind_site_u_F=1

    #4 encode user germline_somatic
    if germline_somatic_u=='germline':
        germline_somatic_u_g=1
        germline_somatic_u_s=0
    else:
        germline_somatic_u_g=0
        germline_somatic_u_s=1
    #5 construct user feature
    user_features=pd.Series([SIFT_score_u,PPH_score_u,GERP_score_u,\
                             volume_change_index_u,polarity_change_index_u,
                             allele_freq_u,\
                             molecular_weight_u,H_bond_donor_u,H_bond_acceptor_u,rotatable_bond_u,Polar_Surface_Area_u,\
                             bind_site_u_F,bind_site_u_T,germline_somatic_u_g,germline_somatic_u_s]).values.reshape(1,-1)
    user_features=pd.DataFrame(user_features)

    training_features_name=['SIFT_score','PPH_score','GERP_score',\
                             'volume_change_index','polarity_change_index',
                             'allele_freq',\
                             'molecular_weight','H_bond_donor','H_bond_acceptor','rotatable_bond','Polar_Surface_Area',\
                             'bind_site_F','bind_site_T','germline_somatic_g','germline_somatic_s']
    return user_features,training_features_name

###################################################################################################################
#                                                   main                                                          #
###################################################################################################################
def main():
    ###################################################################################################################
    #                                                   input_parameter_check                                         #
    ###################################################################################################################
    global mutated_amino_acid_u
    global wildtype_amino_acid_u
    global SIFT_score_u
    global PPH_score_u
    global GERP_score_u
    global allele_freq_u
    global bind_site_u
    global germline_somatic_u
    global molecular_weight_u
    global H_bond_donor_u
    global H_bond_acceptor_u
    global rotatable_bond_u
    global Polar_Surface_Area_u
    ##1 Check mutated_amino_acid
    amino_acid_list=['ALA','LYS','GLU','ASP','ARG','TYR','TRP','THR','SER','HIS',
                     'GLN','CYS','ASN','VAL','PRO','PHE','MET','LEU','ILE','GLY']
    if (mutated_amino_acid_u not in amino_acid_list):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('mutated_amino_acid'))
        sys.exit(1)
    ##2 Check wildtype_amino_acid
    if (wildtype_amino_acid_u not in amino_acid_list):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('wildtype_amino_acid'))
        sys.exit(1)
    ##3 Check SIFT_score
    if SIFT_score_u.replace('.', '', 1).isnumeric():
        SIFT_score_u=float(SIFT_score_u)
        if SIFT_score_u<0 or SIFT_score_u>1:
            print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('SIFT_score'))
            sys.exit(1)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('SIFT_score'))
        sys.exit(1)
    ##4 Check PPH_score
    if PPH_score_u.replace('.', '', 1).isnumeric():
        PPH_score_u=float(PPH_score_u)
        if PPH_score_u<0 or PPH_score_u>1:
            print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit'.format('PPH_score'))
            sys.exit(1)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('PPH_score'))
        sys.exit(1)
    ##5 Check Gerp_score
    if GERP_score_u.replace('.', '', 1).isnumeric():
        GERP_score_u=float(GERP_score_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('GERP_score'))
        sys.exit(1)
    ##6 Check allele_freq
    if allele_freq_u.replace('.', '', 1).isnumeric():
        allele_freq_u=float(allele_freq_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('allele_freq'))
        sys.exit(1)
    ##7 Check bind_site
    if (bind_site_u=='True') or (bind_site_u=='False'):
        bind_site_u=bool(bind_site_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('bind_site'))
        sys.exit(1)
    ##8 Check germline_somatic
    if (germline_somatic_u != 'germline' and germline_somatic_u != 'somatic'):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('germline_somatic'))
        sys.exit(1)
    ##9 Check molecular_weight
    if molecular_weight_u.replace('.', '', 1).isnumeric():
        molecular_weight_u=float(molecular_weight_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('molecular_weight'))
        sys.exit(1)
    ##10 Check H_bond_donor_u
    if H_bond_donor_u.replace('.', '', 1).isnumeric():
        H_bond_donor_u=float(H_bond_donor_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('H_bond_donor'))
        sys.exit(1)
    ##11 Check H_bond_acceptor_u
    if H_bond_acceptor_u.replace('.', '', 1).isnumeric():
        H_bond_acceptor_u=float(H_bond_acceptor_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('H_bond_acceptor'))
        sys.exit(1)
    ##12 Check rotatable_bond_u
    if rotatable_bond_u.replace('.', '', 1).isnumeric():
        rotatable_bond_u=float(rotatable_bond_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('rotatable_bond'))
        sys.exit(1)
    ##13 Check allele_freq
    if Polar_Surface_Area_u.replace('.', '', 1).isnumeric():
        Polar_Surface_Area_u=float(Polar_Surface_Area_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('Polar_Surface_Area'))
        sys.exit(1)
    ###################################################################################################################
    #                                           finish  input_parameter_check                                         #
    ###################################################################################################################

    clf = joblib.load("trained_random_forest.m")
    user_features,training_features_name=parse_user_input(mutated_amino_acid_u,wildtype_amino_acid_u,\
                        SIFT_score_u,PPH_score_u,GERP_score_u,allele_freq_u,\
                        bind_site_u,germline_somatic_u,\
                        molecular_weight_u,H_bond_donor_u,H_bond_acceptor_u,rotatable_bond_u,Polar_Surface_Area_u)

    pred_proba=clf.predict_proba(np.array(user_features).reshape(1,-1))

    print ('Random Forest Features:')
    for name, value in zip(training_features_name,user_features.values.tolist()[0]):
        print (name,':',value)

    print ('')
    print ('Probability that binding_affinity_score_change_class is positive: ',round(float(pred_proba[:,1]),3)*100,'%')
    print ('Probability that binding_affinity_score_change_class is negative or zero: ',round(float(pred_proba[:,0]),3)*100,'%')

main()
