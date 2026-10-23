#!/usr/bin/env Rscript
# Adapted from Ref_Meth_hg38/scripts/rds_to_bed.R.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script <- gsub("~+~", " ", script, fixed = TRUE) # Rscript encodes spaces in --file.
source(file.path(dirname(normalizePath(script)), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: rds_to_bed.R REFERENCE_RDS OUTPUT_DIRECTORY MIN_SAMPLES")
ref <- read_reference(args[1])
out <- args[2]
minimum <- as.integer(args[3])
x <- copy(ref$sites)
pnames <- paste0("p", seq(5, 95, 5))
metrics <- c(pnames, "hits", "coverage_pct", "score", "iqr")
x[, hits := lengths(values)]
if (nrow(x)) {
  qs <- t(vapply(x$values, function(v) quantile(v, probs = seq(.05, .95, .05),
                   na.rm = TRUE, names = FALSE, type = 7), numeric(19)))
  for (i in seq_along(pnames)) set(x, j = pnames[i], value = qs[, i])
  x[, coverage_pct := 100 * hits / length(ref$sample_ids)]
  x[, score := 100 * hits / max(hits)]
  x[, iqr := p75 - p25]
} else {
  for (name in setdiff(metrics, "hits")) set(x, j = name, value = numeric())
}
x[, values := NULL]
setcolorder(x, c("chr", "start", "end", metrics))
dir.create(file.path(out, "statistics"), recursive = TRUE, showWarnings = FALSE)
stats_path <- file.path(out, "statistics", paste0(ref$chromosome, ".tsv.gz"))
if (nrow(x)) {
  fwrite(x, stats_path, sep = "\t", na = "NA")
} else {
  # Explicitly close an empty gzip stream; some fwrite versions omit its trailer.
  con <- gzfile(stats_path, "wt")
  writeLines(paste(names(x), collapse = "\t"), con)
  close(con)
}
reported <- x[hits >= minimum]
for (metric in metrics) {
  dest <- file.path(out, "per_chromosome", metric)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  fwrite(reported[, c("chr", "start", "end", metric), with = FALSE],
         file.path(dest, paste0(ref$chromosome, ".bedGraph")), sep = "\t",
         col.names = FALSE, quote = FALSE, na = "NA")
}
dir.create(file.path(out, "qc"), recursive = TRUE, showWarnings = FALSE)
summary <- data.table(chromosome = ref$chromosome, observed_loci = nrow(x),
                      reported_loci = nrow(reported), cohort_samples = length(ref$sample_ids),
                      median_hits = if (nrow(x)) median(x$hits) else NA_real_)
fwrite(summary, file.path(out, "qc", paste0(ref$chromosome, ".tsv")), sep = "\t", na = "NA")
