
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

library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(Rsamtools)
library(devtools)
library(biovizBase)
library(tidyverse)
library(RCurl)
library(stringr)
library(optparse)


# Command line arguments
option_list = list(
  make_option(
    c("-i", "--in_fragments_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The fragments+gc file to get bin compartments for. Must have the extension .rds."
  ),
  make_option(
    c("-b", "--in_bin_coordinates_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A file with bin coordinates. Required when `--assembly` is 'hg38'. Must have the extension .rds."
  ),
  make_option(
    c("-o", "--out_bins_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The output file containing the 100kb bins. Must have the extension `.rds`."
  ),
  make_option(
    c("-f", "--filters_dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The directory containing the `filters.hg38.rda` and `gaps.hg38.rda` files (or their hg19 counterparts)."
  ),
  make_option(
    c("-s", "--sample_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The ID of the sample."
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
if (is.na(opt$in_fragments_file)){
  stop("--in_fragments_file was not specified.")
}
if (is.na(opt$out_bins_file)){
  stop("--out_bins_file was not specified.")
}
if (is.na(opt$filters_dir)){
  stop("--filters_dir was not specified.")
}
if (is.na(opt$sample_id)){
  stop("--sample_id was not specified.")
}
if (is.na(opt$assembly)){
  stop("--assembly was not specified.")
}
if (str_sub(opt$in_fragments_file, start= -4) != ".rds"){
  stop("--in_fragments_file must have the extension '.rds'.")
}
if (!is.na(opt$in_bin_coordinates_file) && str_sub(opt$in_bin_coordinates_file, start= -4) != ".rds"){
  stop("--in_bin_coordinates_file must have the extension '.rds'.")
}
if (str_sub(opt$out_bins_file, start= -4) != ".rds"){
  stop("--out_bins_file must have the extension '.rds'.")
}
if (!(opt$assembly %in% c("hg19", "hg38"))){
  stop("--assembly was one of {'hg19', 'hg38'}.")
}

print("Got these command line arguments: ")
print(opt)

if (opt$assembly == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
  load(file.path(opt$filters_dir, "filters.hg19.rda")); filters <- filters.hg19
  load(file.path(opt$filters_dir, "gaps.hg19.rda")); gaps <- gaps.hg19

  if (is.na(opt$in_bin_coordinates_file)){
    ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    AB <- read.table(textConnection(ABurl), header = TRUE)
  } else {
    AB <- readRDS(file = opt$in_bin_coordinates_file)
  }
} else {
  library(BSgenome.Hsapiens.UCSC.hg38)
  load(file.path(opt$filters_dir, "filters.hg38.rda")); filters <- filters.hg38
  load(file.path(opt$filters_dir, "gaps.hg38.rda")); gaps <- gaps.hg38
  AB <- readRDS(file = opt$in_bin_coordinates_file)
}

print("Loaded AB file")

AB <- makeGRangesFromDataFrame(
  AB,
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = TRUE # Added this
)

print("Made GRanges object")

gc.correct <- function(coverage, bias) {
  i <- seq(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias)
  coverage.model <- loess(predict(coverage.trend, i) ~ i)
  coverage.pred <- predict(coverage.model, bias)
  coverage.corrected <-
    coverage - coverage.pred + median(coverage)
}

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))

tcmeres <- gaps[grepl("centromere|telomere", gaps$type)]

arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
print("Non-tel/centromere arms"); print(arms)
arms <- arms[-c(25,27,29,41,43)]
print("Filtered arms"); print(arms)

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

arms$arm <- armlevels

print("Set these arm levels: ")
print(arms$arm)

AB <- AB[-queryHits(findOverlaps(AB, gaps))]
print("Removed intervals overlapping with gaps")
AB <- AB[queryHits(findOverlaps(AB, arms))]
print("Extracted intervals overlapping with chromosome arms")

print("AB:")
print(AB)

AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]

seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(Hsapiens, AB)

if (opt$assembly == "hg19") {
  ## These bins had no coverage in hg19
  # TODO Check this for hg38?
  AB <- AB[-c(8780, 13665)]
}
fragments <- readRDS(opt$in_fragments_file)

### Filters
fragments <- fragments[-queryHits(findOverlaps(fragments, filters))]
w.all <- width(fragments)

fragments <- fragments[which(w.all >= 100 & w.all <= 220)]
w <- width(fragments)

frag.list <- split(fragments, w)

counts <- sapply(frag.list, function(x)
  countOverlaps(AB, x))
if (min(w) > 100) {
  m0 <- matrix(
    0,
    ncol = min(w) - 100,
    nrow = nrow(counts),
    dimnames = list(rownames(counts), 100:(min(w) - 1))
  )
  counts <- cbind(m0, counts)
}

olaps <- findOverlaps(fragments, AB)
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))
bingc <- rep(NA, length(bin.list))
bingc[unique(subjectHits(olaps))] <-
  sapply(bin.list, function(x)
    mean(x$gc))

### Get modes
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
modes <- Mode(w)
medians <- median(w)
q25 <- quantile(w, 0.25)
q75 <- quantile(w, 0.75)

short <- rowSums(counts[, 1:51])
long <- rowSums(counts[, 52:121])
ratio <- short / long
short.corrected = gc.correct(short, bingc)
long.corrected = gc.correct(long, bingc)
nfrags.corrected = gc.correct(short + long, bingc)
ratio.corrected = gc.correct(ratio, bingc)

AB$short <- short
AB$long <- long
AB$ratio <- short / long
AB$nfrags <- short + long
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc

for (i in 1:ncol(counts))
  elementMetadata(AB)[, colnames(counts)[i]] <- counts[, i]

# Convert to tibble
AB <- AB %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(id = sample_id)

saveRDS(AB, opt$out_bins_file)
q('no')
