
# Select assembly version
genome <- "hg38" # or hg19

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
  options(timeout = 1000)
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
  options(timeout = 1000) # Might fail due to slow downloads
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
}
if (!requireNamespace("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")


library(GenomicRanges)
if (genome == "hg19") {
  library(BSgenome.Hsapiens.UCSC.hg19)
} else {
  library(BSgenome.Hsapiens.UCSC.hg38)
}
library(rtracklayer)
library(tidyverse)
library(httr)

if (genome == "hg19"){

  ### hg19 gaps & blacklisted regions
  mySession <- browserSession()
  genome(mySession) <- genome
  gaps <- getTable(ucscTableQuery(mySession, table="gap"))
  gaps.hg19 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                       gaps$chromEnd),
                   type=gaps$type)
  gaps.hg19 <- keepSeqlevels(gaps.hg19, paste0("chr", c(1:22, "X", "Y")),
                             pruning.mode="coarse")
  hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  seqinfo(gaps.hg19) <- seqinfo(hsapiens)[seqlevels(gaps.hg19),]
  save(gaps.hg19, file="gaps.hg19.rda")
  # devtools::use_data(gaps.hg19, overwrite = TRUE)

  blacklisted.file <- httr::content(GET("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
  blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
                              col_names=c("seqnames", "start",
                                          "end", "name", "score"))
  blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
  filters.hg19 <- makeGRangesFromDataFrame(blacklisted.tib,
                                             keep.extra.columns=TRUE)
  filters.hg19 <- keepSeqlevels(filters.hg19, paste0("chr", c(1:22, "X", "Y")),
                             pruning.mode="coarse")
  seqinfo(filters.hg19) <- seqinfo(Hsapiens)[seqlevels(filters.hg19),]
  save(filters.hg19, file="filters.hg19.rda")
  # devtools::use_data(filters.hg19, overwrite = TRUE)

} else if (genome == "hg38"){

  ### hg38 gaps & blacklisted regions
  mySession <- browserSession()
  genome(mySession) <- genome
  gaps <- getTable(ucscTableQuery(mySession, table="gap")) %>%
    dplyr::select(chrom, chromStart, chromEnd, type)
  # Centromere coordinates were extracted following https://www.biostars.org/p/435003/#462800
  centromeres_hg38 <- read.csv("centromere_coordinates_hg38.tsv", sep="\t") %>%
    dplyr::mutate(type = "centromere") %>%
    dplyr::select(chrom, chromStart, chromEnd, type) %>%
    dplyr::group_by(chrom, type) %>%
    # Collapse the two coordinate sets per chrom
    dplyr::summarise(
      chromStart = min(chromStart),
      chromEnd = max(chromEnd),
      .groups = "drop"
    )

  gaps <- dplyr::bind_rows(gaps, centromeres_hg38) %>%
    dplyr::arrange(chrom, chromStart)

  gaps.hg38 <- GRanges(
    gaps$chrom,
    IRanges(gaps$chromStart, gaps$chromEnd),
    type=gaps$type
  )
  gaps.hg38 <- keepSeqlevels(
    gaps.hg38, paste0("chr", c(1:22, "X", "Y")),
    pruning.mode="coarse"
  )
  hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  seqinfo(gaps.hg38) <- seqinfo(hsapiens)[seqlevels(gaps.hg38),]
  save(gaps.hg38, file="gaps.hg38.rda")
  # devtools::use_data(gaps.hg38, overwrite = TRUE)

  blacklisted.file <- httr::content(GET("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"))
  blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
                              col_names=c("seqnames", "start",
                                          "end", "name", "score"))
  blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
  filters.hg38 <- makeGRangesFromDataFrame(blacklisted.tib,
                                           keep.extra.columns=TRUE)
  filters.hg38 <- keepSeqlevels(filters.hg38, paste0("chr", c(1:22, "X", "Y")),
                                pruning.mode="coarse")
  seqinfo(filters.hg38) <- seqinfo(Hsapiens)[seqlevels(filters.hg38),]
  save(filters.hg38, file="filters.hg38.rda")
  # devtools::use_data(filters.hg38, overwrite = TRUE)

} else {
  stop("Did not recognize `genome`")
}
