
#
# Extract fragments and GC contents
#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GenomicAlignments", quietly = TRUE))
  BiocManager::install("GenomicAlignments")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")
if (!requireNamespace("Rsamtools", quietly = TRUE))
  BiocManager::install("Rsamtools")
if (!requireNamespace("Homo.sapiens", quietly = TRUE))
  BiocManager::install("Homo.sapiens")
if (!requireNamespace("biovizBase", quietly = TRUE))
  BiocManager::install("biovizBase")
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(devtools)
library(Homo.sapiens)
library(biovizBase)
library(optparse)


# Command line arguments
option_list = list(
  make_option(
    c("-i", "--in_galp_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The GALP file to get fragments and GC contents for. Must have the extension .rds."
  ),
  make_option(
    c("-o", "--out_fragments_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file for fragments and GC contents in. Must have the extension `.rds`."
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
if (is.na(opt$in_galp_file)){
  stop("--in_galp_file was not specified.")
}
if (is.na(opt$out_fragments_file)){
  stop("--out_fragments_file was not specified.")
}
if (is.na(opt$assembly)){
  stop("--assembly was not specified.")
}
if (str_sub(opt$in_galp_file, start= -4) != ".rds"){
  stop("--in_galp_file must have the extension '.rds'.")
}
if (str_sub(opt$out_fragments_file, start= -4) != ".rds"){
  stop("--out_fragments_file must have the extension '.rds'.")
}
if (!(opt$assembly %in% c("hg19", "hg38"))){
  stop("--assembly was one of {'hg19', 'hg38'}.")
}

if (opt$assembly == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
} else {
  library(BSgenome.Hsapiens.UCSC.hg38)
}

galp <- readRDS(opt$in_galp_file)
frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode = "coarse"),
                 on.discordant.seqnames = "drop")

## Filter outliers
w.all <- width(frags)
q.all <- quantile(w.all, c(0.001, 0.999))
frags <- frags[which(w.all > q.all[1] & w.all < q.all[2])]

gcs <- GCcontent(Hsapiens, unstrand(frags))
frags$gc <- gcs

saveRDS(frags, opt$out_fragments_file)
q('no')
