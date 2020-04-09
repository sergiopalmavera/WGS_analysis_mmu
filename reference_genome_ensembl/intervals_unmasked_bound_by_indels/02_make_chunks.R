library(dplyr)
library(GenomicRanges)

# Import sorted intervals (unmasked-indels-subtracted & indels)
intervals <- read.table("reference_genome_ensembl/intervals_unmasked_bound_by_indels/intervals_sorted.bed",stringsAsFactors = F)
head(intervals)
dim(intervals) #19362081        3

# Calculute rolling sum
intervals <- intervals %>% 
  dplyr::rename(chrom = V1, start = V2, end = V3) %>% 
  group_by(chrom) %>% #group by chrom
  # calculate cummulative lenth of intervals
  mutate(size = end - start, rolling_sum = cumsum(size)) 

head(intervals)
dim(intervals) #19362081        5

# Define maximun size of interval
max_size <- 5000000

# Iterate over each chromosome (only autosomes and X)
intervals_bounded_by_indels <- lapply(c(1:19,"X"), function(chr){
  #chr=3
  # Subset data by chromosome
  d <- intervals %>% filter(chrom == chr)
  # Empty vector to store chunk levels
  levs <- vector(mode = "character", length = nrow(d))
  # Make desired ranges of max_size at a time
  col1 <- seq(from = 0, to = max(d$rolling_sum), by = max_size)
  col2 <- seq(from = max_size, to = max(d$rolling_sum), by = max_size) %>% c(max(.)+max_size)
  # Store desired ranges in a df
  dd <- data.frame(start = col1+1, end = col2) #add a 1 to start to avoid overlaps
  # Populate chunk levels
  for(i in 1:nrow(dd)){
    # classify each entry in rolling_sum according to its chunk/range
    # Each interval is classified depending of the range is contained in
    levs[d$rolling_sum >= dd$start[i] & d$rolling_sum <= dd$end[i]] <- paste0("chunk_",i)
  }
  # split rolling_sum into chunks spanning no more than max_size
  chunks <- split(d$rolling_sum, levs)
  # Iterate over each chunk 
  res <- lapply(unique(levs), function(chunk_i){
    # subset data by rolling_sum values in chunk
    d_chunk_i <- d[d$rolling_sum %in% chunks[[chunk_i]],]
    # Convert to a granges
    gr_i <- GRanges(seqnames = d_chunk_i$chrom, 
                    IRanges(start = d_chunk_i$start, end = d_chunk_i$end))
    # Reduce granges (if gap >1L, don't merge)
    gr_i_red <- reduce(gr_i, min.gapwidth=1L)
    # make list a single granges
   }) %>% purrr::reduce(c)
})

intervals_bounded_by_indels <- intervals_bounded_by_indels %>% purrr::reduce(c)
intervals_bounded_by_indels %>% length() #929

# Export 
intervals_bounded_by_indels %>% 
  as.data.frame() %>% 
  dplyr::select(seqnames, start, end) %>% 
  write.table("./reference_genome_ensembl/intervals_unmasked_bound_by_indels/intervals_bounded_by_indels_max5Mb.bed", col.names = F, row.names = F, quote = F)

