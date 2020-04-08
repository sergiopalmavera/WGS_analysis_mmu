library(dplyr)
library(GenomicRanges)

# Import indel intervals
indels <- read.table("reference_genome_ensembl/intervals_indels/mus_musculus_indels_sorted_lite.bed", 
                     stringsAsFactors = F) %>% filter(V1 %in% c(1:19,"X"))
head(indels)
dim(indels) #10443650        3
# Import unmasked intervals after indels were removed
unmasked_excl_indels <- read.table("reference_genome_ensembl/unmasked_intervals/intervals_unmasked_excl_indels.bed", 
                                   stringsAsFactors = F) %>% filter(V1 %in% c(1:19,"X"))
head(unmasked_excl_indels)
dim(unmasked_excl_indels) #8956927       3
# Combine unmasked-excl indels with indels into a single intervals set
intervals <- bind_rows(unmasked_excl_indels, indels) %>%
  arrange(V1,V2) %>% # sort by chr and start of interval
  dplyr::rename(chrom = V1, start = V2, end = V3) %>% 
  group_by(chrom) %>% #group by chrom
  # calculate cummulative lenth of intervals
  mutate(size = end - start, rolling_sum = cumsum(size)) 

head(intervals)
dim(intervals) #19400379        3
# Export sorted intervals
write.table(intervals, 
            "./reference_genome_ensembl/intervals_unmasked_bound_by_indels/intervals.bed",
            col.names = F, row.names = F, quote = F)
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
    # Each interval is classified depending of the range is contained
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

# Range of intervals
intervals_bounded_by_indels %>% width() %>% sort() %>% range()

# Export 
intervals_bounded_by_indels %>% 
  as.data.frame() %>% 
  dplyr::select(seqnames, start, end) %>% 
  write.table("./reference_genome_ensembl/intervals_unmasked_bound_by_indels/intervals_bounded_by_indels_max5Mb.bed", 
              col.names = F, row.names = F, quote = F)

# Sanity check: are all rows in intervals preserved? --> YES!
nrow(intervals) #19400379

sanity <- lapply(c(1:19,"X"), function(chr){
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
    # Each interval is classified depending of the range is contained
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
    #gr_i_red <- reduce(gr_i, min.gapwidth=1L)
    # make list a single granges
  }) %>% purrr::reduce(c)
})

sanity %>% purrr::reduce(c) %>% length() #19400379

# Are the last rows preserved? #--> yes

xx <- sanity %>% lapply(function(x) x[length(x)]) %>% 
  purrr::reduce(c) %>% as.data.frame() %>% 
  mutate(id = paste(seqnames,start,end,sep = "_")) %>% 
  dplyr::select(id)

yy <- lapply(c(1:19,"X"), function(chr){
  d <- intervals %>% filter(chrom == chr)
  d[nrow(d),] %>% 
    mutate(id = paste(chrom,start,end,sep = "_")) %>% 
    dplyr::select(id)
}) %>% bind_rows()

identical(xx$id, yy$id)




