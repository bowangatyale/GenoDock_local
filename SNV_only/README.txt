test dir: /rf_without_distance


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