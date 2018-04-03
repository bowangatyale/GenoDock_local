test dir: /rf_without_distance


Parameter input in order:
#1 mutated_amino_acid,2 wildtype_amino_acid
#3 SIFT_score,4 PPH_score,5 GERP_score,6 allele_freq,7 bind_site
#8 germline_somatic,9 molecular_weight,10 H_bond_donor
#11 H_bond_acceptor,12 rotatable_bond,13 Polar_Surface_Area

Example1:
1
Input:
python train_random_forest.py
Output:
Start Train Random Forest Model.
Finish Train Random Forest Model.

Get trained model file trained_random_forest.m

2
Input:
python test.py THR GLY 0.05 0.016 4.09 0.0000165 False germline 272.09 5 8 4 153

Output:
Random Forest Features:
SIFT_score : 0.05
PPH_score : 0.016
GERP_score : 4.09
volume_change_index : 0.9541963103868752
polarity_change_index : 0.5
allele_freq : 1.65e-05
molecular_weight : 272.09
H_bond_donor : 5.0
H_bond_acceptor : 8.0
rotatable_bond : 4.0
Polar_Surface_Area : 153.0
bind_site_F : 1.0
bind_site_T : 0.0
germline_somatic_g : 1.0
germline_somatic_s : 0.0

Probability that binding_affinity_score_change_class is positive:  0.2 %
Probability that binding_affinity_score_change_class is negative or zero:  99.8 %

