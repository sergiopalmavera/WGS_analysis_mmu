# Intro
This directory contains all the steps to replicate the analysis of the DU mouse lines WGS data 

# The original data: FASTQ files in 3 sets.

- Set1: the original batch, 60 samples (10 for each line), 9 lanes per sample, paired, 60 * 9 * 2 = 1080 FASTQ files

- Set2: complementary set to boost low coverage samples (n=10), 1 lane per sample, paired, 10 * 2 = 20 FASTQ files

- Set3: new samples to increase N (10/line -> 25/line). 90 samples (15 per line), 1 lane per sample, paired = 90 * 2 = 180 FASTQ files

# The repository is organized in four main parts:

- Part1: "set1" and "set2" data processing from FASTQ to (post-GenotypeRefinement)VCF ==> VCF1

- Part2: "set3" data processing from FASTQ to (post-GenotypeRefinement)VCF ==> VCF2

- Part3: VCF1 and VCF2 combined and filtered to produce the final VCF with 25 * 6 = 150 samples. ==> VCFfinal

- Part4: Population Structure Analysis 

- Part5: Selection Analysis.

# Misc

- There was a drop-out sample (I34772). It was resequenced. The original drop-out was not included in the analysis, the re-sequenced sample was included instead.

- Note that scripts that were run in chinks (i.e. samples 1 to 10) correspond to the order of samples before the drop-out reseqed was added. Keep this in mind and modify scripts accordingly.
