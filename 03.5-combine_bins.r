
library(tidyverse)
library(GenomicRanges)
library(optparse)
library(stringr)

# Command line arguments
option_list = list(
  make_option(
    c("-i", "--in_bin_files"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The comma-separated paths to the bin files. Each path must have the extension `.rds`"
  ),
  make_option(
    c("-o", "--out_bins_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file containing the concatenated 100kb bins. Must have the extension `.rds`."
  )
)
opt = parse_args(OptionParser(option_list = option_list))
if (is.na(opt$in_bin_files)){
  stop("--in_bin_files was not specified.")
}
if (is.na(opt$out_bins_file)){
  stop("--out_bins_file was not specified.")
}
# Extract paths
bin_files <- unlist(strsplit(opt$in_bin_files, ","))
for (bf in bin_files){
  if (str_sub(bf, start= -4) != ".rds"){
    stop("paths in --in_bin_files must have the extension '.rds'.")
  }
}
if (str_sub(opt$out_bins_file, start= -4) != ".rds"){
  stop("--out_bins_file must have the extension '.rds'.")
}

all_bins <- lapply(bin_files, readRDS) %>%
  bind_rows() %>%
  select(id, everything()) %>%
  select(-matches("X"))

saveRDS(all_bins, opt$out_bins_file)
q('no')
