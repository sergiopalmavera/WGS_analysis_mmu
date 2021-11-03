# Intro
This directory contains all the steps to replicate the analysis of the DU mouse lines WGS data. 

# The original data: FASTQ files in 3 sets.
- Batch1: the original batch, 60 samples (10 for each line), 9 lanes per sample, paired, 60 * 9 * 2 = 1080 FASTQ files. Target coverage 30x.
- Batch2: complementary set to boost low coverage samples (n=10), 1 lane per sample, paired, 10 * 2 = 20 FASTQ files. Target coverage 30x.
- Batch3: new samples to increase N (10/line -> 25/line). 90 samples (15 per line), 1 lane per sample, paired = 90 * 2 = 180 FASTQ files. Target coverage 5x.

# Directories
- "00_dashboard": it contains the R markdown scripts to produce dashboards to visualize and summarize the main parts of the analysis.
- "01_scripts_batch1_and_batch2": Scripts to process the reads (01_fastQC), alignments (02_bwa), prepare bams and run haplotype-caller (03_bamprep_bqsr_haplotypecaller) on each sample in batch1 and batch2. It also contains scripts used to collect metrics from alignments.
- "batch3" same as "01_scripts_batch1_and_batch2" but for batch3 samples.
- "batches123_01_ConsolidateGVFs_by_intervals" scripts to consolidate gVCFs produced by haplotype caller.
- "batches123_02_GenotypeGVCFs_by_chr" estimate genotypes from consolidated VCF.
- "batches123_03_VariantQualityScoreRecalibration" scripts to recalibrate quality scores at each variant site.
- "batches123_04_final_VCF" filtering of raw VCF and metric collection.
- "batches123_05_annotation" annotation of variants with SNPeff and postprocessing of results.
- "batches123_06_genetic_structure": PCA, hierarchical clustering and ADMIXTURE analysis.
- "batches123_07_selection_analysis": Fst and pi calculations. Not the best name for directory. 
- "batches123_08_LDD_SFS": analysis of linkage disequilibrium decay based on genotipic correlations between SNPs, and allele frequency visualizations.
- "batches123_09_ROH": estimate runs of homozygosity with bcftools roh.
- "reference_genome_ensembl": reference genome and associated files downloaded from ensembl 93.
- "resource_mgp_sanger_REL-1505-SNPs_Indels": processing of MGP variant data used as resource for variant quality score recalibration.
- "sample_info" sample information
- "selection_experiment" scripts to process genotypes of control animals to select them according to their genetic similarity to the fertility lines.

# Data accessibility:
- Raw reads at the European Nucleotide Archive (accession: PRJEB44248)
- VCFs at the European Variation Archive (accession: PRJEB45961)

# Manuscript
- The preprint for this study can be found at: https://www.biorxiv.org/content/10.1101/2021.05.28.446207v1.full
- A simplified version of these repository: https://github.com/sosfert/mmu_dummerstorf_wgs

# Misc
- There was a drop-out sample (I34772). It was resequenced. The original drop-out was not included in the analysis, the re-sequenced sample was included instead.
- Original drop out sample was excluded from downstream analysis. Only resequenced sample was used.
