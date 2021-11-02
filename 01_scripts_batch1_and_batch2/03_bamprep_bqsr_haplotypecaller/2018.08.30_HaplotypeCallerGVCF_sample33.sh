# Print hostname
echo "Hostname"
hostname

# Define path raw fastqs
RAWFASTQ=/projekte/I2-SOS-FERT/Original

# Set path to directory with BQSR-BAMs
BAM_BQSR=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/results/merged

# Define path to GATK
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

# Define path to temprary files
TMP=/projekte/I2-SOS-FERT/05_HaplotypeCaller_GVCF/tmp

# Define path for output
OUT=/projekte/I2-SOS-FERT/05_HaplotypeCaller_GVCF/results

# Define sample to process using original fastq directory as reference (consistent with former scripts)
FROM=33
TO=33

SAMPLES=$(ls -1 $RAWFASTQ | sed '/md5/d' | grep .fastq.gz | cut -d'_' -f1 | sort | uniq | sed -n $FROM,${TO}p)

time for sample in $SAMPLES
do
	echo "## SAMPLE: $sample"
	time samtools index $BAM_BQSR/$sample.merged.sorted.dedup.bqsr.bam 
	mv $BAM_BQSR/$sample.merged.sorted.dedup.bqsr.bam.bai $BAM_BQSR/$sample.merged.sorted.dedup.bqsr.bai
	time $GATK/gatk HaplotypeCaller \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $BAM_BQSR/$sample.merged.sorted.dedup.bqsr.bam \
		-O $TMP/$sample.g.vcf.gz \
		-ERC GVCF
	mv $TMP/$sample.g.vcf.gz $OUT
done

