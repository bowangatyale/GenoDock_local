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

#dir
PDB_FILE_DIR = "PDB/"
PDB_CONFIG_FILE = "manually_annotated_BLASTClump_output.csv"

BIND_SITE_OUTPUT_DIR = "BindSite/"
RESIDUE_DISTANCE_OUTPUT_DIR = "ResidueDistance/"
RESIDUE_DISTANCE_DIR = "ResidueDistance/"
#hide warning
pd.options.mode.chained_assignment = None  # default='warn'

###################################################################################################################
#                                              input parameter read                                               #
###################################################################################################################
#1 read command line argument
#Parameters users input:
#1 pdb_id,2 chain_id,3 ligand_u,4 amino_acid_position,5 mutated_amino_acid,
#6 wildtype_amino_acid,7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
#11 germline_somatic
argument_list=sys.argv
pdb_id_u=argument_list[1]
chain_id_u=argument_list[2]
ligand_u=argument_list[3]
amino_acid_position_u=argument_list[4]
mutated_amino_acid_u=argument_list[5]
wildtype_amino_acid_u=argument_list[6]
SIFT_score_u=argument_list[7]
PPH_score_u=argument_list[8]
GERP_score_u=argument_list[9]
allele_freq_u=argument_list[10]
germline_somatic_u=argument_list[11]
molecular_weight_u=argument_list[12]
H_bond_donor_u=argument_list[13]
H_bond_acceptor_u=argument_list[14]
rotatable_bond_u=argument_list[15]
Polar_Surface_Area_u=argument_list[16]



###################################################################################################################
#                                              compute residue distance                                           #
###################################################################################################################
DISTANCE_BOUNDARY = 8

def sort_atom_set_key(atom):
    return int(atom[1])

def sort_aa_set_key(aa):
    return int(aa[1])
def parse_pdb_config(pdb_config_filename):
    config = {}
    with open(pdb_config_filename) as csv_file:
        reader = csv.reader(csv_file, delimiter=',', quotechar='|')
        next(reader)
        for row in reader:
            config[row[0]] = (row[1], row[3])

    return config
def calculate_residue_distance(pdb_config):
    min_distance = 999999

    for file in os.listdir(PDB_FILE_DIR):

        file_path = PDB_FILE_DIR + file
        filename = os.path.splitext(file)[0]
        ext = os.path.splitext(file)[1]

        if ext == ".pdb" and filename==pdb_id_u:

            search_ligand_id = pdb_config[file][1]
            search_chain = pdb_config[file][0]

            #print("Parsing {}, Chain: {}, Ligand: {}......".format(filename, search_chain, search_ligand_id))

            pdb_file = open(file_path, "r")
            result_file = open(RESIDUE_DISTANCE_OUTPUT_DIR + filename + ".txt", "w")

            atom_coordinates = []
            hatatm_coordinates = []
            amino_acid_set = set()
            amino_acid_distance_min = {}

            ligand_id_exist=[]
            for pdbLineInfo in pdb_file:

                if pdbLineInfo[0:6] == 'ATOM  ':

                    chain = pdbLineInfo[21].strip()
                    if chain != search_chain:
                        continue

                    x = pdbLineInfo[31:38].strip()
                    y = pdbLineInfo[40:46].strip()
                    z = pdbLineInfo[48:54].strip()

                    amino_acid_name = pdbLineInfo[16:20].strip()[-3:]
                    position = pdbLineInfo[22:26].strip()
                    amino_acid_distance_min[(amino_acid_name, position, chain)] = min_distance

                    atom_coordinates.append((x, y, z, amino_acid_name, position, chain))

                if pdbLineInfo[0:6] == 'HETATM':
                    ligand_id = pdbLineInfo[17:21].strip()
                    if ligand_id == search_ligand_id:
                        ligand_id_exist.append(ligand_id)
                        x = pdbLineInfo[31:38].strip()
                        y = pdbLineInfo[40:46].strip()
                        z = pdbLineInfo[48:54].strip()
                        hatatm_coordinates.append((x, y, z))

            if search_ligand_id not in ligand_id_exist:
                return 0

            for atom in atom_coordinates:
                atom_x = float(atom[0])
                atom_y = float(atom[1])
                atom_z = float(atom[2])

                key = (atom[3], atom[4], atom[5])

                for hetatm in hatatm_coordinates:
                    hetatm_x = float(hetatm[0])
                    hetatm_y = float(hetatm[1])
                    hetatm_z = float(hetatm[2])

                    x_square = (atom_x - hetatm_x) * (atom_x - hetatm_x)
                    y_square = (atom_y - hetatm_y) * (atom_y - hetatm_y)
                    z_square = (atom_z - hetatm_z) * (atom_z - hetatm_z)
                    distance = (x_square + y_square + z_square) ** (1. / 2.)


                    if distance < amino_acid_distance_min[key]:
                        amino_acid_distance_min[key] = distance

            for k in amino_acid_distance_min:
                distance = amino_acid_distance_min[k]
                if distance != min_distance:
                    amino_acid_set.add((k[0], k[1], k[2], distance))

            for aa in sorted(amino_acid_set, key=sort_aa_set_key):
                result_file.write(aa[0] + "\t" + aa[1] + '\t' + aa[2] + '\t' + str(aa[3]) + "\n")

            result_file.close()
            return 1
###################################################################################################################
#                                                   parse_user_features                                                  #
###################################################################################################################
#Randomforest model required features:
#1 SIFT_score,2 PPH_score,3 GERP_score,4 volumn_change_index
#5 distance,6 distance_normed,7 allele_freq,
#8 bind_site_F,9 bind_site_T,10 germline_somatic_g,11 germline_somatic_s

#User Input
#1 pdb_id,2 chain_id,3 ligand_u,4 amino_acid_position,5 mutated_amino_acid,
#6 wildtype_amino_acid,7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
#11 germline_somatic 12 molecular_weight,13 H_bond_donor,14 H_bond_acceptor,
#15rotatable_bond,16 Polar_Surface_Area
def parse_user_input(pdb_id_u,chain_id_u,amino_acid_position_u,\
                    mutated_amino_acid_u,wildtype_amino_acid_u,\
                    SIFT_score_u,PPH_score_u,GERP_score_u,allele_freq_u,\
                    germline_somatic_u,molecular_weight_u,H_bond_donor_u,\
                    H_bond_acceptor_u,rotatable_bond_u,Polar_Surface_Area_u):
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

    #3 compute distance normed
    from sklearn.preprocessing import MinMaxScaler
    PDB_FILE_DIR = "PDB/"
    RESIDUE_DISTANCE_DIR = "ResidueDistance/"

    #####construct a distance dictionary
    distance={}
    for pdb_file in os.listdir(RESIDUE_DISTANCE_DIR):
        protein_name=pdb_file[0:4]
        if pdb_file[0] !='.' and protein_name==pdb_id_u:
            protein_name=pdb_file[0:4]
            p_path=os.path.join('ResidueDistance',pdb_file)
            pdistance=pd.read_table(p_path,index_col=False,header=None)[3].values.reshape(-1,1)
            scaler=MinMaxScaler()
            distance[protein_name]=pd.Series(list(scaler.fit_transform(pdistance).transpose()[0]))
    ####concat protein txt files with 2 more columns: pdb_id column, normed distance column
    i=0
    for pdb_file in os.listdir("ResidueDistance/"):
        #print (pdb_file)
        protein_name=pdb_file[0:4]
        if pdb_file[0] !='.' and protein_name==pdb_id_u:
            #####read every distance txt file
            p_path=os.path.join('ResidueDistance',pdb_file)

            df=pd.read_table(p_path,index_col=False,header=None)
            #####add pdb_id column
            row_num=df.shape[0]
            df['pdb_id']=row_num*[protein_name]
            #####add normed distance column
            df['distance_normed']=pd.Series(distance[protein_name])
    #####concat all df
    i=i+1
    if i==1:
        distance_concat=df
    else:
        distance_concat=pd.concat([distance_concat,df])

    distance_concat.columns=['wildtype_amino_acid','amino_acid_position',
                            'chain_id','distance','pdb_id','distance_normed']
    distance_concat=distance_concat.drop('wildtype_amino_acid',1)
    #4 compute bindsite according to distance: distance<8,bs=True
    bind_site_u=[]
    for distance in distance_concat['distance']:
        if distance<8:
            bind_site_u.append(True)
        else:
            bind_site_u.append(False)
    distance_concat['bind_site']=bind_site_u
    distance_bindsite_concat=distance_concat.reset_index(drop=True)
    distance_bindsite_concat.to_csv('distance_bindsite_concat.csv')

    #5 construct user feature
    user_df=pd.Series([pdb_id_u,chain_id_u,amino_acid_position_u,mutated_amino_acid_u,\
                         wildtype_amino_acid_u,SIFT_score_u,PPH_score_u,GERP_score_u,allele_freq_u,\
                         volume_change_index_u,polarity_change_index_u,germline_somatic_u,molecular_weight_u,H_bond_donor_u,\
                         H_bond_acceptor_u,rotatable_bond_u,Polar_Surface_Area_u]).values.reshape(1,-1)
    user_df=pd.DataFrame(user_df,
                columns=['pdb_id','chain_id','amino_acid_position','mutated_amino_acid',\
                         'wildtype_amino_acid','SIFT_score','PPH_score','gerp_score','allele_freq',\
                         'volume_change_index','polarity_change_index','germline_somatic','molecular_weight','H_bond_donor',\
                         'H_bond_acceptor','rotatable_bond','Polar_Surface_Area'])
    ###merge to add distance/distance_normed/bind_site
    user_df=user_df.merge(distance_bindsite_concat,\
                            how='inner',\
                            on=['pdb_id','chain_id',\
                                'amino_acid_position'])

    #6 encode user bind_site
    bind_site_T=[]
    bind_site_F=[]
    for value in user_df['bind_site']:
        if value==True:
            bind_site_T.append(1)
            bind_site_F.append(0)
        else:
            bind_site_T.append(0)
            bind_site_F.append(1)
    user_df['bind_site_F']=bind_site_F
    user_df['bind_site_T']=bind_site_T


    # encode user germline_somatic
    germline_somatic_g=[]
    germline_somatic_s=[]
    for value in user_df['germline_somatic']:
        if value=='germline':
            germline_somatic_g.append(1)
            germline_somatic_s.append(0)
        else:
            germline_somatic_g.append(0)
            germline_somatic_s.append(1)
    user_df['germline_somatic_g']=germline_somatic_g
    user_df['germline_somatic_s']=germline_somatic_s

    user_features=user_df[['SIFT_score', 'PPH_score', 'gerp_score', 'volume_change_index','polarity_change_index',\
                          'distance', 'distance_normed', 'allele_freq', 'molecular_weight', 'H_bond_donor', \
                          'H_bond_acceptor', 'rotatable_bond', 'Polar_Surface_Area', 'bind_site_F', 'bind_site_T', \
                          'germline_somatic_g', 'germline_somatic_s']]
    return user_features.loc[0],list(user_features.columns)

###################################################################################################################
#                                                   main                                                          #
###################################################################################################################
def main():
    ###################################################################################################################
    #                                                   input_parameter_check                                         #
    ###################################################################################################################
    global pdb_id_u
    global chain_id_u
    global ligand_u
    global amino_acid_position_u
    global mutated_amino_acid_u
    global wildtype_amino_acid_u
    global SIFT_score_u
    global PPH_score_u
    global GERP_score_u
    global allele_freq_u
    global germline_somatic_u
    global molecular_weight_u
    global H_bond_donor_u
    global H_bond_acceptor_u
    global rotatable_bond_u
    global Polar_Surface_Area_u

    ##1 Check pdb_id
    if (len(pdb_id_u)!=4):
        print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('pdb_id'))
        sys.exit(1)
    ##2 Check chain_id
    if (len(chain_id_u)!=1 or chain_id_u.isupper()=='False'):
        print ('Input parameter {} is invaliad, pleasd double check the input format and resubmit.'.format('chain_id'))
        sys.exit(1)
    ##3 check ligand id
    #   at line 479

    ##4 Check amino_acid_position
    if amino_acid_position_u.isnumeric():
        amino_acid_position_u=int(amino_acid_position_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('amino_acid_position'))
        sys.exit(1)
    ##5 Check mutated_amino_acid
    amino_acid_list=['ALA','LYS','GLU','ASP','ARG','TYR','TRP','THR','SER','HIS',
                     'GLN','CYS','ASN','VAL','PRO','PHE','MET','LEU','ILE','GLY']
    if (mutated_amino_acid_u not in amino_acid_list):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('mutated_amino_acid'))
        sys.exit(1)
    ##6 Check wildtype_amino_acid
    if (wildtype_amino_acid_u not in amino_acid_list):
        print ('Input {} is invaliad, pleasd double check the input format and resubmit.'.format('wildtype_amino_acid'))
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
    ##12 Check molecular_weight
    if molecular_weight_u.replace('.', '', 1).isnumeric():
        molecular_weight_u=float(molecular_weight_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('molecular_weight'))
        sys.exit(1)

    ##13 Check H_bond_donor_u
    if H_bond_donor_u.replace('.', '', 1).isnumeric():
        H_bond_donor_u=float(H_bond_donor_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('H_bond_donor'))
        sys.exit(1)
    ##14 Check H_bond_acceptor_u
    if H_bond_acceptor_u.replace('.', '', 1).isnumeric():
        H_bond_acceptor_u=float(H_bond_acceptor_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('H_bond_acceptor'))
        sys.exit(1)
    ##15 Check rotatable_bond_u
    if rotatable_bond_u.replace('.', '', 1).isnumeric():
        rotatable_bond_u=float(rotatable_bond_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('rotatable_bond'))
        sys.exit(1)
    ##16 Check allele_freq
    if Polar_Surface_Area_u.replace('.', '', 1).isnumeric():
        Polar_Surface_Area_u=float(Polar_Surface_Area_u)
    else:
        print ('Input {} is invaliad, pleasd double check the input format and resubmit'.format('Polar_Surface_Area'))
        sys.exit(1)
    ###################################################################################################################
    #                                           finish  input_parameter_check                                         #
    ###################################################################################################################


    #1compute distance
    config=np.array([pdb_id_u+'.pdb',chain_id_u,amino_acid_position_u,ligand_u]).reshape(1,-1)
    input_config=pd.DataFrame(config,columns=['PDB','mapped_chain','amino_acid_position','ligand'])
    input_config.to_csv('manually_annotated_BLASTClump_output.csv',index=False)
    pdb_config = parse_pdb_config(PDB_CONFIG_FILE)
    calculate_result=calculate_residue_distance(pdb_config)

    if calculate_result==0:
        print ('Input {} is invaliad, cannot find ligand_id in your provided pdb file'.format('ligand_id'))
        sys.exit(1)


    #load rf model
    clf = joblib.load("trained_random_forest.m")
    #3parse_user_input
    user_features,training_features_name=parse_user_input(pdb_id_u,chain_id_u,amino_acid_position_u,\
                                                          mutated_amino_acid_u,wildtype_amino_acid_u,\
                                                          SIFT_score_u,PPH_score_u,GERP_score_u,allele_freq_u,\
                                                          germline_somatic_u,\
                                                          molecular_weight_u,H_bond_donor_u,H_bond_acceptor_u,rotatable_bond_u,Polar_Surface_Area_u)
    #4predict
    pred_proba=clf.predict_proba(np.array(user_features).reshape(1,-1))
    #5print result
    print ('Random Forest Features:')
    for name, value in zip(training_features_name,user_features.values.tolist()):
        print (name,':',value)

    print ('')
    print ('Probability that binding_affinity_score_change_class is positive: ',round(float(pred_proba[:,1]),3)*100,'%')
    print ('Probability that binding_affinity_score_change_class is negative or zero: ',round(float(pred_proba[:,0]),3)*100,'%')
main()
