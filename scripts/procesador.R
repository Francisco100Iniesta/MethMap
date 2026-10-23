#!/usr/bin/env Rscript
# Adapted from the original chromosome splitter. Paths and contigs are now inputs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: procesador.R NORMALIZED_SAMPLE_DIRECTORY")
suppressPackageStartupMessages(library(data.table))
setDTthreads(as.integer(Sys.getenv("METHMAP_R_THREADS", "1")))
files <- list.files(args[1], pattern = "\\.tsv$", full.names = TRUE)
for (path in files) {
  x <- fread(path, sep = "\t", colClasses = c(chr = "character"))
  if (!identical(names(x), c("chr", "start", "end", "percent", "coverage"))) stop("Invalid input columns")
  if (anyDuplicated(x[, .(chr, start, end)])) {
    stop("Duplicate coordinates in ", path, "; supply one unphased 5mC record per cytosine and sample.")
  }
  setorder(x, start, end)
  saveRDS(x, sub("\\.tsv$", ".rds", path), compress = "gzip")
  unlink(path)
}
