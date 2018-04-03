Command Line Input Parameters:
#1 pdb_id,2 chain_id,3 ligand_u,4 amino_acid_position,5 mutated_amino_acid,
#6 wildtype_amino_acid,7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
#11 germline_somatic 12 molecular_weight,13 H_bond_donor,14 H_bond_acceptor,
#15rotatable_bond,16 Polar_Surface_Area


Example:
1
Input:
python train_random_forest.py
Output:
Start Train Random Forest Model.
Finish Train Random Forest Model.

Get trained model file trained_random_forest.m

2
Input:
python test.py 1t46 A STI 687 ASN GLY 0.2 0.601 4.67 0.00004 germline 272.09 5 8 4 153

Output:
Random Forest Features:
SIFT_score : 0.2
PPH_score : 0.601
gerp_score : 4.67
volume_change_index : 1.0
polarity_change_index : 0.5
distance : 20.4128153129
distance_normed : 0.623218185803
allele_freq : 4e-05
molecular_weight : 272.09
H_bond_donor : 5.0
H_bond_acceptor : 8.0
rotatable_bond : 4.0
Polar_Surface_Area : 153.0
bind_site_F : 1
bind_site_T : 0
germline_somatic_g : 1
germline_somatic_s : 0

Probability that binding_affinity_score_change_class is positive:  0.1 %
Probability that binding_affinity_score_change_class is negative or zero:  99.9 %
