
# Added based on https://github.com/cancer-genomics/reproduce_lucas_wflow/issues/1#issuecomment-1055679129
# NOTE: Requires a LOT of RAM. So set R_MAX_VSIZE=100Gb or similar in the environment. 16gb was not enough.

if (!requireNamespace("optparse", quietly = TRUE))
  BiocManager::install("optparse")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
if (!requireNamespace("Homo.sapiens", quietly = TRUE))
  BiocManager::install("Homo.sapiens")

library(optparse)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Homo.sapiens)

# Command line arguments
option_list = list(
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output `cytosine_ref_{assembly}.rds` file for the specified assembly."
  ),
  make_option(
    c("-a", "--assembly"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The genomic assembly to use. Either 'hg19' or 'hg38'."
  )
)
opt = parse_args(OptionParser(option_list = option_list))

if (is.na(opt$out_file)){
  stop("--out_file was not specified.")
}
if (is.na(opt$assembly)){
  stop("--assembly was not specified.")
}
if (str_sub(opt$out_file, start= -4) != ".rds"){
  stop("--out_file must have the extension '.rds'.")
}
if (!(opt$assembly %in% c("hg19", "hg38"))){
  stop("--assembly was one of {'hg19', 'hg38'}.")
}

if (opt$assembly == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
} else {
  library(BSgenome.Hsapiens.UCSC.hg38)
}

cytosines <- vmatchPattern("C", Hsapiens)

saveRDS(cytosines, opt$out_file)
