
library(optparse)
library(tidyverse)
library(purrr)

# Command line arguments
option_list = list(
  make_option(
    c("-i", "--chrom_sizes_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A file mapping the chromosome name to its size. Will be read with `read_tsv` and should have no header."
  ),
  make_option(
    c("-o", "--out_bins_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file containing coordinates for the 100kb bins. Must have the extension `.rds`."
  ),
  make_option(
    c("-s", "--bin_size"),
    action = "store",
    default = 100000,
    type = 'integer',
    help = "The bin size."
  )
)
opt = parse_args(OptionParser(option_list = option_list))
if (is.na(opt$chrom_sizes_file)){
  stop("--chrom_sizes_file was not specified.")
}
if (is.na(opt$out_bins_file)){
  stop("--out_bins_file was not specified.")
}

chrom_sizes_df <- read_tsv(opt$chrom_sizes_file, col_names = c("chr", "size")) %>%
  dplyr::filter(chr %in% paste0("chr", 1:22))
chrom_sizes <- chrom_sizes_df[["size"]]
names(chrom_sizes) <- chrom_sizes_df[["chr"]]

create_chromsome_bins <- function(chrom_name, end, start=0, bin_size=100000){
  num_bins <- floor((end-start)/bin_size)
  start_coordinates <- (seq_len(num_bins) - 1) * bin_size + start
  end_coordinates <- start_coordinates + bin_size - 1
  tibble::tibble(
    "chr" = chrom_name,
    "start" = start_coordinates,
    "end" = end_coordinates
  )
}

create_bin_coordinates <- function(chrom_sizes, bin_size = 100000) {
  purrr::map2_df(.x = chrom_sizes,
                 .y = names(chrom_sizes),
                 .f = ~ {
                   create_chromsome_bins(chrom_name = .y,
                                         end = .x,
                                         bin_size = bin_size)
                 }) %>%
    dplyr::arrange(chr)
}

bin_coordinates <- create_bin_coordinates(
  chrom_sizes=chrom_sizes,
  bin_size=opt$bin_size
)

saveRDS(bin_coordinates, file=opt$out_bins_file)
q('no')
