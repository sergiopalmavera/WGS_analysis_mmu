library(tidyverse)
library(here)
library(vroom)
library(GenomicRanges)
library(stringr)
library(matrixStats)

# Prelude
### Genotype matrix
FZTDU_Fwd_SNP_matrix <- here("/selection experiment/SNP-selection/","FZTDU_Fwd_SNP-matrix.txt") %>% vroom()
FZTDU_Fwd_SNP_matrix

### Marker information
Marker_Giga_MUGA <- here("/selection experiment/SNP-selection","Marker_Giga_MUGA.csv") %>% vroom()
Marker_Giga_MUGA %>% nrow()


#However, there are some markers with strange chromosomes (0?, NA?)
Marker_Giga_MUGA %>% group_by(Chr) %>% summarise(n()) %>% as.data.frame()

### Cross mapping build 37 to build 38
#GIGAMUGA markers ("Marker_Giga_MUGA.csv") are based on build version 37, they were converted to coordinates in build version 38 (ENSEMBL: "assembly converter"). This is necessary for compatibility with WGS analysis

load("snps.gigamuga.Rdata")  #http://csbio.unc.edu/MUGA/
freq_snps<-read.csv("SNP_Freq_FZTDU.csv")
freq_snps<-subset(freq_snps,freq_snps$Minor.Freq>=0.05)

snps<-merge(freq_snps[1],snps[,1:6],by.x="Name",by.y="marker",sort=F)
snps$chr <- sub("chr","",snps$chr)

Marker_Giga_MUGA_buid_37_to_38<-snps
Marker_Giga_MUGA_buid_37_to_38_GR <- Marker_Giga_MUGA_buid_37_to_38 %>% 
filter(!is.na(chr)) %>% 
mutate(tmp = pos) %>% 
dplyr::rename(seqnames = chr, start = pos, end = tmp) %>% 
dplyr::select(seqnames, start, end, everything()) %>% 
makeGRangesFromDataFrame(keep.extra.columns = TRUE)

animal_info <- here("/selection experiment/SNP-selection/","DNA_FZTDU_Gen_37_SNP_chip2.csv") %>% 
vroom() %>% 
mutate(animal_id = paste0("mmu_",`lfd. Nr.`)) %>% 
dplyr::select(animal_id, everything())

#animal_info
animal_info %>% group_by(Sex) %>% summarise(n())

# Subset genotypes according to regions of distinct genetic differentiation

win_fst_sel_vs_ctrl <- lapply(
list.files(here("Z:/FBNGroup/I4-SOS-FERT/03_SAW/LinuX/projekte_I2-SOS-FERT/GitHub/WGS_analysis_mmu/batches123_07_selection_analysis/01_win_fst_sel_vs_ctrl/output"),
pattern = "windowed.weir.fst$"), 

function(fl){
read.table(
here("Z:/FBNGroup/I4-SOS-FERT/03_SAW/LinuX/projekte_I2-SOS-FERT/GitHub/WGS_analysis_mmu/batches123_07_selection_analysis/01_win_fst_sel_vs_ctrl/output",fl),
header = T,
stringsAsFactors = F
)
}
)

names(win_fst_sel_vs_ctrl) <- list.files(
here("Z:/FBNGroup/I4-SOS-FERT/03_SAW/LinuX/projekte_I2-SOS-FERT/GitHub/WGS_analysis_mmu/batches123_07_selection_analysis/01_win_fst_sel_vs_ctrl/output"),
pattern = "windowed.weir.fst$"
) %>% str_remove(".windowed.weir.fst")



# prepare data (i.e. zscore transformation) for visualizations
win_fst_sel_vs_ctrl_prep <- lapply(win_fst_sel_vs_ctrl, function(d){
d %>% 
# filter by min number of SNPs in window (min 10 SNPs)
dplyr::filter(N_VARIANTS >= 10) %>% 
# separate x chromosome
mutate(is_x = CHROM == "X") %>% 
# group by autosomes/chrX
group_by(is_x) %>% 
# compute zscores
mutate(z_win_fst_score = scale(MEAN_FST))   
}) %>% 
bind_rows(.id = "contrast") %>% 
mutate(contrast = factor(
contrast, 
levels = c("DUK_vs_FZTDU","DUC_vs_FZTDU","DU6_vs_FZTDU","DU6P_vs_FZTDU","DUhLB_vs_FZTDU","FERT_vs_FZTDU"))
)

# Import quatniles
summary_z_win_fst <- readRDS( here("Z:/FBNGroup/I4-SOS-FERT/03_SAW/LinuX/projekte_I2-SOS-FERT/GitHub/WGS_analysis_mmu/00_dashboard","summary_z_win_fst.rds") ) # produced elsewhere (check 00_dashboard.Rmd)

### Define functions to extract highest and lowest Fst quantiles
get_top_windows <- function(contr, qtile){
qtlie_zscore <- summary_z_win_fst[summary_z_win_fst$contrast == contr, qtile] %>% as.numeric()
win_fst_sel_vs_ctrl_prep %>% 
ungroup() %>% 
filter(contrast == contr) %>% 
filter(z_win_fst_score > qtlie_zscore) %>% 
dplyr::select(CHROM, BIN_START, BIN_END) 
}

get_bottom_windows <- function(contr, qtile){
qtlie_zscore <- summary_z_win_fst[summary_z_win_fst$contrast == contr, qtile] %>% as.numeric()
win_fst_sel_vs_ctrl_prep %>% 
ungroup() %>% 
filter(contrast == contr) %>% 
filter(z_win_fst_score < qtlie_zscore) %>% 
dplyr::select(CHROM, BIN_START, BIN_END) 
}

### Extract regions of distinct genetic differentiation (also refered to as signatures by prof Reinsch)

chrs_ordered <- c(1:19,"X")

RDD_DUK_q15_vs_q85 <- get_top_windows("DUK_vs_FZTDU", "q85") %>% 
inner_join(get_bottom_windows("DUC_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
inner_join(get_bottom_windows("DU6_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
inner_join(get_bottom_windows("DU6P_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
inner_join(get_bottom_windows("DUhLB_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
# merge adjacent and overlapping windows
reduce(min.gapwidth = 2000000) %>% 
# define the name of the signature
as_tibble() %>% 
mutate(signature_id = paste0(seqnames,":",start, "-", end)) %>% 
# return to granges
makeGRangesFromDataFrame(keep.extra.columns = TRUE)

RDD_DUC_q15_vs_q85 <- get_top_windows("DUC_vs_FZTDU", "q85") %>% 
inner_join(get_bottom_windows("DUK_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
inner_join(get_bottom_windows("DU6_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
inner_join(get_bottom_windows("DU6P_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
inner_join(get_bottom_windows("DUhLB_vs_FZTDU","q15"), by = c("CHROM","BIN_START","BIN_END")) %>% 
mutate(CHROM = factor(CHROM, levels = chrs_ordered)) %>% 
mutate(CHROM = factor(CHROM, levels = chrs_ordered),
signature_id = paste0(CHROM,":",BIN_START, "-", BIN_END)) %>% 
makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
# merge adjacent and overlapping windows
reduce(min.gapwidth = 2000000) %>% 
# define the name of the signature
as_tibble() %>% 
mutate(signature_id = paste0(seqnames,":",start, "-", end)) %>% 
# return to granges
makeGRangesFromDataFrame(keep.extra.columns = TRUE)


prep_overlaps <- function(gr_gigamuga, gr_sig){
# find overlaps
ov <- findOverlaps(query = gr_gigamuga, subject = gr_sig, maxgap = 1000000)
# extract overlaps from gigamuga
gr_gigamuga_ov <- gr_gigamuga[queryHits(ov)]

# extract overalaps from RDDs
gr_sig_ov <- gr_sig[subjectHits(ov)]

# add signature information to gigamuga overlaps
res <- bind_cols(
gr_gigamuga_ov %>% as_tibble() %>% dplyr::select(-width, -strand),
gr_sig_ov %>% as_tibble() %>% dplyr::select(signature_id)
)
return(res)
}

gigamuga_x_RDDq15q85_DUK <- prep_overlaps(Marker_Giga_MUGA_buid_37_to_38_GR, RDD_DUK_q15_vs_q85)
#write.csv(gigamuga_x_RDDq15q85_DUK,"DUK_sel.csv",row.names = F)

gigamuga_x_RDDq15q85_DUC <- prep_overlaps(Marker_Giga_MUGA_buid_37_to_38_GR, RDD_DUC_q15_vs_q85)
#write.csv(gigamuga_x_RDDq15q85_DUC,"DUC_sel.csv",row.names = F)

gigamuga_x_RDDq15q85_DUK %>% nrow()
gigamuga_x_RDDq15q85_DUC %>% nrow()
gigamuga_x_RDDq15q85_DUK$signature_id %>% unique() %>% length()
gigamuga_x_RDDq15q85_DUC$signature_id %>% unique() %>% length()

# Convert genotypes to allele counts
# get allele freq
af_prep <- function(pop){
  path_nm <- paste0("Z:/FBNGroup/I4-SOS-FERT/03_SAW/LinuX/projekte_I2-SOS-FERT/GitHub/WGS_analysis_mmu/batches123_08_LDD_SFS/output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.allrecords.",pop,".ALTfrq2_fixed_cols")
  path_nm %>% vroom(col_names = FALSE, col_types = c("X1"="c")) %>% dplyr::rename(CHROM=X1, POS=X2)
} 

af_dat <- vroom(
  file="Z:/FBNGroup/I4-SOS-FERT/03_SAW/LinuX/projekte_I2-SOS-FERT/GitHub/WGS_analysis_mmu/batches123_04_final_VCF/output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.SNPids2", 
  col_names = FALSE, 
  delim = " ",
  col_types = c(X1 = "c")
)%>% 
  dplyr::select(-X2,-X4,-X6) %>% 
  dplyr::rename(CHROM = X1, POS = X3, REF=X5,ALT=X7) %>% 
  inner_join(af_prep("FZTDU"), by = c("CHROM", "POS")) %>% 
  inner_join(af_prep("DUK"), by = c("CHROM", "POS")) %>% 
  inner_join(af_prep("DUC"), by = c("CHROM", "POS"))
lines_sorted<-c("FZTDU","DUK","DUC")
names(af_dat)[grep("^X",names(af_dat))] <- lines_sorted


overlap_af <- inner_join(af_dat, snps, by = c("CHROM"="chr","POS"="pos"))# %>%   #overlap between frequent (FZTDU) SNPs and Chip-SNPs
  mutate(ref_match = (REF == A1))

#geno <- paste0("SNP-selection/FZTDU_Fwd_SNP-matrix.txt") %>% vroom()
#marker<- Marker_Giga_MUGA_buid_37_to_38_GR[,c("marker","A1","A2")]
#mg<-merge(marker,geno,by.x="marker",by.y="SNP Name",sort=F)
#save(mg,file="markergeno.RDa")
load("markergeno.RDa")
rownames(mg)<-mg[,1]

#for DUK
geno.DUK<-merge(gigamuga_x_RDDq15q85_DUK,mg,by.x="Name",by.y=0)
geno.DUC<-merge(gigamuga_x_RDDq15q85_DUC,mg,by.x="Name",by.y=0)

geno.DUK.af<-merge(overlap_af[,c(8,5:7)],geno.DUK, by="Name",sort=F)
geno.DUK.af.Ref<-subset(geno.DUK.af,geno.DUK.af$DUK<=geno.DUK.af$FZTDU)
geno.DUK.af.Alt<-subset(geno.DUK.af,geno.DUK.af$DUK>=geno.DUK.af$FZTDU)

num.DUK.Ref<- apply(geno.DUK.af.Ref[,grep("^1$",colnames(geno.DUK.af.Ref)):dim(geno.DUK.af.Ref)[2]],2,str_count,pattern=geno.DUK.af.Ref$A1.x)
num.DUK.Ref<-cbind(geno.DUK.af.Ref[,c(1,2,3,8)],num.DUK.Ref)

num.DUK.Alt<- apply(geno.DUK.af.Alt[,grep("^1$",colnames(geno.DUK.af.Alt)):dim(geno.DUK.af.Alt)[2]],2,str_count,pattern=geno.DUK.af.Alt$A2.x)
num.DUK.Alt<-cbind(geno.DUK.af.Alt[,c(1,2,3,8)],num.DUK.Alt)
num.DUK<-rbind(num.DUK.Ref,num.DUK.Alt)

DUK<-merge(geno.DUK[,c(1,2,3,8)],num.DUK,by="Name",sort=F)
dim(DUK)

write.csv(DUK,"DUK_sel_SNPs.csv",row.names = F)
DUK.sum<-aggregate(DUK[,grep("^1$",colnames(DUK)):dim(DUK)[2]],by=list(DUK$signature_id),"sum")
DUK.var<-as.data.frame(rowMeans(DUK.sum[,-1]))
colnames(DUK.var)<-"mean"
DUK.var$SD<-rowSds(as.matrix(DUK.sum[,-1]))
DUK.norm<-(DUK.sum[,-1]-DUK.var[,1])/DUK.var[,2]
DUK.norm.anno<-cbind(DUK.sum[,1],DUK.norm)
DUK.CI<-colSums(DUK.norm)
write.csv(DUK.CI,"DUK_mouse_selection.csv")
hist(DUK.CI)
head(DUK.CI[order(DUK.CI)]);tail(DUK.CI[order(DUK.CI)])

#for DUC
geno.DUC.af<-merge(overlap_af[,c(8,5:7)],geno.DUC, by="Name",sort=F)
geno.DUC.af.Ref<-subset(geno.DUC.af,geno.DUC.af$DUC<=geno.DUC.af$FZTDU)
geno.DUC.af.Alt<-subset(geno.DUC.af,geno.DUC.af$DUC>=geno.DUC.af$FZTDU)

num.DUC.Ref<- apply(geno.DUC.af.Ref[,grep("^1$",colnames(geno.DUC.af.Ref)):dim(geno.DUC.af.Ref)[2]],2,str_count,pattern=geno.DUC.af.Ref$A1.x)
num.DUC.Ref<-cbind(geno.DUC.af.Ref[,c(1,2,4,8)],num.DUC.Ref)

num.DUC.Alt<- apply(geno.DUC.af.Alt[,grep("^1$",colnames(geno.DUC.af.Alt)):dim(geno.DUC.af.Alt)[2]],2,str_count,pattern=geno.DUC.af.Alt$A2.x)
num.DUC.Alt<-cbind(geno.DUC.af.Alt[,c(1,2,4,8)],num.DUC.Alt)
num.DUC<-rbind(num.DUC.Ref,num.DUC.Alt)

DUC<-merge(geno.DUC[,c(1,2,3,8)],num.DUC,by="Name",sort=F)
dim(DUC)

write.csv(DUC,"DUC_sel_SNPs.csv",row.names = F)
DUC.sum<-aggregate(DUC[,grep("^1$",colnames(DUC)):dim(DUC)[2]],by=list(DUC$signature_id),"sum")
DUC.var<-as.data.frame(rowMeans(DUC.sum[,-1]))
colnames(DUC.var)<-"mean"
DUC.var$SD<-rowSds(as.matrix(DUC.sum[,-1]))
DUC.norm<-(DUC.sum[,-1]-DUC.var[,1])/DUC.var[,2]
DUC.norm.anno<-cbind(DUC.sum[,1],DUC.norm)
DUC.CI<-colSums(DUC.norm)
write.csv(DUC.CI,"DUC_mouse_selection.csv")
hist(DUC.CI)
head(DUC.CI[order(DUC.CI)]);tail(DUC.CI[order(DUC.CI)])
