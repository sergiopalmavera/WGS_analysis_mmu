#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=$1
max_size=5000000
max_nr_jobs=$2

echo "# Arguments passed:"
echo "## chr = $chr"
echo "## max_nr_jobs = $max_nr_jobs"
printf "\n"

echo "# Running Intervals for Chromosome $chr"

mkdir ../output/chr$chr #Make a chromosome file
mkdir ./chr${chr}_logs # Make a subdir to store log files
mkdir ./intervals/chr$chr # Make a subdirectory for intervals
dir_iv=./intervals/chr$chr # store directory for intervals for simplicity

BED=../../reference_genome_ensembl/unmasked_intervals # path to unmasked intervals
BED_INDELS=../../reference_genome_ensembl/intervals_indels
BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2 #absolute path, modify accordingly

# Get intervals equal or larger than max_size
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked_excl_indels.bed | awk '{print $0,  $3 - $2}' | awk -v ms=$max_size -F' ' '$4 >= ms {print $1 "\t" $2 "\t" $3}' > $dir_iv/chr${chr}_large_i.bed

# Export small intervals
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked_excl_indels.bed | awk '{print $0,  $3 - $2}' | awk -v ms=$max_size -F' ' '$4 < ms {print $1 "\t" $2 "\t" $3}' > $dir_iv/chr${chr}_small_i.bed

# Split large intervals into intervals of max_size
$BEDTOOLS/bedtools.static.binary makewindows -b $dir_iv/chr${chr}_large_i.bed -w $max_size > $dir_iv/chr${chr}_large_i_split.bed

# Concatenate small intervals with large intervals split
cat $dir_iv/chr${chr}_small_i.bed $dir_iv/chr${chr}_large_i_split.bed > $dir_iv/intervals_chr${chr}.bed

# Count number of intervals after indels were removed
n_intervals=$(wc -l $dir_iv/intervals_chr${chr}.bed | awk '{print $1}')
echo "## Total nr of intervals excl indels: $n_intervals"

# Split intervals into list of max_nr_jobs
n_lines=$(( (n_intervals / max_nr_jobs) + 1 )) #add another line to make sure no interval is excluded
echo "## Number of intervals per job: $n_lines (total jobs running in parallel = $max_nr_jobs)"

# split intervals into max_nr_jobs chunks of n_lines each
split -l $n_lines $dir_iv/intervals_chr${chr}.bed $dir_iv/chunk --additional-suffix .bed --numeric-suffixes=1

# sanity check
echo "## Sanity check: Number of intervals across jobs (should be same as nr of intervals = $n_intervals)"
cat $dir_iv/chunk* | wc -l
printf "\n"

echo "## Number of batches (also number of GenomicsDBImport instances):"
ls -1 $dir_iv/chunk* | wc -l
printf "\n"

# loop over each interval of max_size (excl indels) for chromosome 
for batch in $dir_iv/chunk*
do
	echo "## Consolidate GVCFs for intervals in batch $batch"	
	nm=$(basename $batch | sed 's/.bed//g')
	nohup ./02_ConsolidateGVCFs_template.sh $batch ../output/chr$chr/$nm &> ./chr${chr}_logs/ConsolidateGVCFs_${nm}.out & # submit job
	printf "\n"
done

# Get indels for chromosome and submit GenomicsDBImport for those intervals
echo "## Consolidate indels as one job (not in parallel)"

awk -v x=$chr '$1 == x {print $1 "\t" $2 "\t" $3}' $BED_INDELS/mus_musculus_indels_sorted.bed > $dir_iv/mus_musculus_indels_sorted_chr${chr}.bed

nohup ./02_ConsolidateGVCFs_template.sh $dir_iv/mus_musculus_indels_sorted_chr${chr}.bed ../output/chr$chr/indels &> ./chr${chr}_logs/ConsolidateGVCFs_${chr}_indels.out &

printf "\n"

echo "## All jobs for chromosome submitted"
printf "\n"


