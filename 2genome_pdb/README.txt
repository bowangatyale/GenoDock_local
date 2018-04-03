 Test dir: /rf_with_distance


Command Line Input Parameters:
#1 pdb_id,2 chain_id,3 amino_acid_position,
#4 mutated_amino_acid,5 wildtype_amino_acid,6 bind_site
#7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
#11 germline_somatic

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
python test.py 2c6n A 8 ASN GLY False 0.92 0.661 9.67 0.004 somatic

Output:
Random Forest Features:
SIFT_score : 0.92
PPH_score : 0.661
GERP_score : 9.67
volume_change_index : 1.0
polarity_change_index : 0.5
allele_freq : 0.004
bind_site_F : 1.0
bind_site_T : 0.0
germline_somatic_g : 0.0
germline_somatic_s : 1.0

Probability that binding_affinity_score_change_class is positive:  5.4 %
Probability that binding_affinity_score_change_class is negative or zero:  94.6 %