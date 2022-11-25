
library(tidyverse)
library(multidplyr)
library(GenomicRanges)
library(readxl)
library(optparse)


# Command line arguments
option_list = list(
  make_option(
    c("-i", "--in_bins_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The input file containing the concatenated 100kb bins for all samples. Each path must have the extension `.rds`"
  ),
  make_option(
    c("-o", "--out_bins_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file containing the summarized 5mb bins. Must have the extension `.rds`."
  ),
  make_option(
    c("-o", "--out_features_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file containing the 5mb bin features. Must have the extension `.csv`."
  )
)
opt = parse_args(OptionParser(option_list = option_list))
if (is.na(opt$in_bins_file)){
  stop("--in_bins_file was not specified.")
}
if (is.na(opt$out_bins_file)){
  stop("--out_bins_file was not specified.")
}
if (str_sub(opt$in_bins_file, start= -4) != ".rds"){
  stop("--in_bins_file must have the extension '.rds'.")
}
if (str_sub(opt$out_bins_file, start= -4) != ".rds"){
  stop("--out_bins_file must have the extension '.rds'.")
}


df.fr <- readRDS(opt$in_bins_file)
#master <- read_csv("sample_reference.csv")

#df.fr2 <- inner_join(df.fr, master, by=c("id"="WGS ID"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr$arm <- factor(df.fr$arm, levels = armlevels)

## combine adjacent 100kb bins to form 5mb bins. We count starting from
## the telomeric end and remove the bin closest to the centromere if it is
## smaller than 5mb.
df.fr2 <- df.fr %>%
  group_by(id, arm) %>%
  mutate(combine = ifelse(
    grepl("p", arm),
    ceiling((1:length(arm)) / 50),
    ceiling(rev((1:length(arm)) / 50))
  ))

# Keep: nfrags.corrected2, short.corrected2, sample, bin,

df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
  summarize(
    #short2=sum(short),
    #long2=sum(long),
    short.corrected2 = sum(short.corrected),
    #long.corrected2=sum(long.corrected),
    #hic.eigen=mean(eigen),
    #gc=mean(C.G),
    #ratio2=mean(ratio),
    #ratio.corrected2=mean(ratio.corrected),
    #nfrags2=sum(nfrags),
    nfrags.corrected2 = sum(nfrags.corrected),
    #domain = median(as.integer(domain)),
    #short.var=var(short.corrected),
    #long.var=var(long.corrected),
    #nfrags.var=var(nfrags.corrected),
    #mode_size=unique(mode),
    #mean_size=unique(mean),
    #median_size=unique(median),
    #q25_size=unique(quantile.25),
    #q75_size=unique(quantile.75),
    #start=start[1],
    #end=rev(end)[1],
    binsize = n()
  ) %>%
  filter(binsize == 50) %>%
  # NOTE: Here they used "sample" instead of id, but that seems like a mistake
  # (don't know where 'sample' would be created and 'id' refers to sample id)
  group_by(id) %>%
  mutate(bin = 1:length(id))

# Save 5mb bins
saveRDS(df.fr3, opt$out_bins_file) # "bins_5mbcompartments.rds"

# Convert to features
features.cov <- df.fr3 %>%
  ungroup() %>%
  select(nfrags.corrected2, id, bin) %>%
  spread(id, nfrags.corrected2) %>%
  select(-bin) %>%
  na.omit() %>%
  scale() %>%
  t() %>%
  as.data.frame()

features.short <- df.fr3 %>%
  ungroup() %>%
  select(short.corrected2, id, bin) %>%
  spread(sample, short.corrected2) %>%
  select(-bin) %>%
  na.omit() %>%
  scale() %>%
  t() %>%
  as.data.frame()

features.sl <- cbind(features.cov, features.short)
colnames(features.sl) <-
  c(paste0("total", 1:498), paste0("short", 1:498))

# Save features in csv
write.csv(features.sl, file = opt$out_features_file)
q('no')
