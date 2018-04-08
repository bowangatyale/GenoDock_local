# User Manual for Genodock Program Suite


## Introduction
GenoDock is a rapid and efficient program suite that could prioritize nsSNVs that would potentially distablize drug ligand-protein binding activity. The software package consists of 4 models based on various feature availability from user. Using SNV annotation, drug ligand molecule, and protein structure features, Genodock will feedback the probability of target nsSNVs to positivly shift bidning affinity between protein and drug ligand.  Here we provide GenoDOck source script for users. For each model, we procide an user manual below.  

### SNV Only
Parameter input in order:
1SIFT_score，2PPH_score,3GERP_score,
4bind_site,5allele_freq,6germline_somatic

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
python test.py 0.05 0.016 4.09 False 0.0000165 germline

Output:
Random Forest Features:
SIFT_score : 0.05
PPH_score : 0.016
GERP_score : 4.09
allele_freq : 1.65e-05
bind_site_F : 1.0
bind_site_T : 0.0
germline_somatic_g : 1.0
germline_somatic_s : 0.0

Probability that binding_affinity_score_change_class is positive:  0.3 %
Probability that binding_affinity_score_change_class is negative or zero:  99.7 %

### SNV + Ligand
Parameter input in order:
1 mutated_amino_acid,2 wildtype_amino_acid
3 SIFT_score,4 PPH_score,5 GERP_score,6 allele_freq,7 bind_site
8 germline_somatic,9 molecular_weight,10 H_bond_donor
11 H_bond_acceptor,12 rotatable_bond,13 Polar_Surface_Area

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

### SNV + PDB
Command Line Input Parameters:
1 pdb_id,2 chain_id,3 amino_acid_position,
4 mutated_amino_acid,5 wildtype_amino_acid,6 bind_site
7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
11 germline_somatic

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

### SNV + Ligand + PDB
Command Line Input Parameters:
1 pdb_id,2 chain_id,3 ligand_u,4 amino_acid_position,5 mutated_amino_acid,
6 wildtype_amino_acid,7 SIFT_score,8 PPH_score,9 GERP_score,10 allele_freq
11 germline_somatic 12 molecular_weight,13 H_bond_donor,14 H_bond_acceptor,
15rotatable_bond,16 Polar_Surface_Area


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
