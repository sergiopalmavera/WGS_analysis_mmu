library(dplyr)
library(ggplot2)
library(stringr)

fls <- list.files("../output",pattern = ".ld")
fls

dat <- lapply(fls, function(fl){
	        d <- read.table(file.path("../output",fl), header = TRUE, stringsAsFactors = FALSE) 
		d <- mutate(d, dist = BP_B - BP_A)
	        return(d)
		})

names(dat) <- str_remove(fls, "cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered_thinned.recode.") %>% str_remove(".plink.ld")
names(dat)

dat2 <- dat %>% bind_rows(.id = "line")
head(dat2)

dat2$line <- factor(dat2$line, levels = c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU"))

# Full data set (max separation any 2 SNPs = 5Mb)
ggplot(data = dat2, aes(x=dist/1000, y=R2, color = line)) +
  geom_smooth(method = "auto") +
    xlab("Pairwise distance (Kb)") +
    ylab(expression(LD~(r^{2}))) +
    theme_bw()+
    theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5)) +
    ggtitle("Linkage Disequilibrium Decay (5Mb)") +
    ggsave("../figures/trends_max5Mb.png")

# Zoom in to 250K
dat2_ss <- dat2 %>% filter(dist <= 250000)

ggplot(data = dat2_ss, aes(x=dist/1000, y=R2, color = line)) +
  geom_smooth(method = "auto") +
  xlab("Pairwise distance (Kb)") +
  ylab(expression(LD~(r^{2}))) +
  theme_bw()+
  theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Linkage Disequilibrium Decay (250Kb)") +
  ggsave("../figures/trends_max250Kb.png")
