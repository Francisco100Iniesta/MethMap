#!/usr/bin/env Rscript
# Adapted from Ref_Meth_hg38/scripts/Add_info.R; writes a new version, never the base.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script <- gsub("~+~", " ", script, fixed = TRUE) # Rscript encodes spaces in --file.
source(file.path(dirname(normalizePath(script)), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Add_info.R BASE_RDS MANIFEST OUTPUT_RDS")
ref <- read_reference(args[1])
samples <- read_manifest(args[2])
for (i in seq_len(nrow(samples))) ref <- append_sample(ref, samples$path[i], samples$sample_id[i])
save_reference(ref, args[3])
