

#
# Extract Genomic Alignment Pairs
#


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GenomicAlignments", quietly = TRUE))
  install.packages("GenomicAlignments")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  install.packages("GenomicRanges")
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")

library(GenomicAlignments)
library(GenomicRanges)
library(optparse)
library(stringr)


# Command line arguments
option_list = list(
  make_option(
    c("-i", "--in_bam_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The bam file to get genomic alignment pairs for."
  ),
  make_option(
    c("-o", "--out_galp_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file for storing the genomic alignment pairs in. Should end in `.rds`."
  ),
)
opt = parse_args(OptionParser(option_list = option_list))
if (is.na(opt$in_bam_file)){
  stop("--in_bam_file was not specified.")
}
if (is.na(opt$out_galp_file)){
  stop("--out_galp_file was not specified.")
}
if (str_sub(opt$out_galp_file, start= -4) != ".rds"){
  stop("--out_galp_file must have the extension '.rds'.")
}

### Read GAlignmentPairs
indexed.bam <- gsub("$", ".bai", opt$in_bam_file)
if (!file.exists(indexed.bam)) {
  indexBam(opt$in_bam_file)
}

param <- ScanBamParam(
  flag = scanBamFlag(
    isDuplicate = FALSE,
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE
  ),
  mapqFilter = 30
)

galp <- readGAlignmentPairs(opt$in_bam_file, param = param)
saveRDS(galp, opt$out_galp_file)
