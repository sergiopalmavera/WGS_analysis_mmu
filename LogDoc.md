mportant #
#############

- I started (quite late I must say) this log book to write down the train of thought I took to solve certain issues. 
- I used to do this on my physical notebook by pen, but since the corona break started Ive been trying to avoid touching objects as much as possible. 
- This document is a dump of mostly errors and how I solve them. Also it is sort of a diary to list the stuff I did on any given day. I guess this might be helpfull for the future.


##################
# Mon 16.03.2020 # 
##################


# GenotypeGVCFs threw an error while processing chr2:

10:08:46.698 INFO  ProgressMeter -          2:153971859           5725.0            1354740000         236636.3
10:08:56.783 INFO  ProgressMeter -          2:154007859           5725.2            1354776000         236635.7
10:09:00.920 INFO  GenotypeGVCFs - Shutting down engine
[15. M▒rz 2020 10:09:01 MEZ] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 5,725.25 minutes.
Runtime.totalMemory()=5758779392
htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 1354789544: there are 120 genotypes while the header requires that 149 genotypes be present for all records at 2:154021236

# At 2:154021236 there are 120 genotypes instead of 149 (150 minus the drop out sample).

# What should I do next?

- Extract the region around the problematic line (this takes a while).

- Post a question on the gatk forum

- After dropout sample is done being processed by HaplotypeCaller and while I wait for an forum-answer, I will use "ConsolidateGVCFs" instead of "CombineGVCFs" (https://gatk.broadinstitute.org/hc/en-us/articles/360035889971?flash_digest=cab2d85e4c16398547fbdae6825867894ff59cf7)

# By the way:
- The number of sites in the cohort.g.vcf file (after CombineGVCFs) is 200,823,577 ... 200M variants!!

- The total number of short variants in dbSNPs is 83,761,978 (http://www.ensembl.org/Mus_musculus/Info/Annotation)

- The size of the file is 4T

- The dropout sample was resequenced and is currently under HaplotypeCaller. Either way I was not gonna move on from GenotypeGVCFs had it been succesfull.




