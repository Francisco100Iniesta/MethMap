#!/usr/bin/env Rscript
# Adapted from Ref_Meth_hg38/scripts/baseconstructor.R: any number of samples.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script <- gsub("~+~", " ", script, fixed = TRUE) # Rscript encodes spaces in --file.
source(file.path(dirname(normalizePath(script)), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: baseconstructor.R CHROMOSOME MANIFEST OUTPUT_RDS")
samples <- read_manifest(args[2])
ref <- list(schema_version = 1L, chromosome = args[1], sample_ids = character(), sites = empty_sites())
for (i in seq_len(nrow(samples))) ref <- append_sample(ref, samples$path[i], samples$sample_id[i])
save_reference(ref, args[3])
