admixture=/home/fb4/palma-vera/FBN_HOME/Tools/ADMIXTURE/admixture_linux-1.3.0

FL=../data/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.ldpruned.bed

$admixture/admixture $FL 2
$admixture/admixture $FL 3
$admixture/admixture $FL 4
$admixture/admixture $FL 5
