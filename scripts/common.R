# MethMap: shared functions adapted from Ref_Meth_hg38 (2026).
suppressPackageStartupMessages(library(data.table))
setDTthreads(as.integer(Sys.getenv("METHMAP_R_THREADS", "1")))

empty_sites <- function() {
  data.table(chr = character(), start = integer(), end = integer(), values = list())
}

save_reference <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Cannot save reference: ", path)
}

read_reference <- function(path) {
  x <- readRDS(path)
  if (!is.list(x) || !identical(x$schema_version, 1L) || is.null(x$sites) ||
      is.null(x$sample_ids) || anyDuplicated(x$sample_ids)) {
    stop("Unsupported reference format. Use a reference built by this MethMap version.")
  }
  x$sites <- as.data.table(x$sites)
  x
}

append_sample <- function(ref, sample_path, sample_id) {
  if (sample_id %in% ref$sample_ids) stop("Duplicate sample ID: ", sample_id)
  x <- readRDS(sample_path)
  if (nrow(x) && any(x$chr != ref$chromosome)) stop("Chromosome mismatch")
  if (anyDuplicated(x[, .(chr, start, end)])) stop("Duplicate sample coordinates")
  x <- x[, .(chr, start, end, percent)]
  if (nrow(x)) {
    joined <- merge(ref$sites, x, by = c("chr", "start", "end"), all = TRUE, sort = FALSE)
    joined$values <- Map(function(old, value) {
      if (is.null(old) || !length(old)) old <- numeric()
      if (is.na(value)) old else c(old, setNames(value, sample_id))
    }, joined$values, joined$percent)
    joined[, percent := NULL]
    setorder(joined, start, end)
    ref$sites <- joined
  }
  ref$sample_ids <- c(ref$sample_ids, sample_id)
  ref
}

read_manifest <- function(path) {
  x <- fread(path, sep = "\t", colClasses = "character")
  if (!identical(names(x), c("sample_id", "path")) || !nrow(x) || anyDuplicated(x$sample_id)) {
    stop("Expected a non-empty manifest with sample_id and path columns")
  }
  x
}
