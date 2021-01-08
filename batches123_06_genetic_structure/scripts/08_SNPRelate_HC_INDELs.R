# Checking example before running SNPRelate on my data
## According to https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate
## To clarify axes labeling: https://support.bioconductor.org/p/119389/

# Load packages
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

genofile <- seqOpen("../data/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.ldpruned.gds")

ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile))

sample_info <- read.csv("../../sample_info/sample_info.csv")

# Make sure that samples have the same order before assigning groups
data.frame(x=ibs.hc$sample.id,y=sample_info$sample_id)

rv <- snpgdsCutTree(ibs.hc, samp.group=as.factor(sample_info$Linie))
rv

saveRDS(rv$dendrogram, "../data/dendrogram_INDELs.rds")

saveRDS(as.factor(sample_info$Linie), "../data/Linie_INDELs.rds")

sessionInfo()
