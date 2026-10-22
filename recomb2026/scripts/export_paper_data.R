# Export every per-feature summary of a results root as one table, for the paper's data release.
#
# Usage: Rscript scripts/export_paper_data.R <out_dir> <dataset> <results_root> [<wald_comparison_dir>] [--methods=<regex>]
#   --methods: export only the method directories whose name matches the regex (e.g. '^nbsr' or '^(?!nbsr)', Perl
#              syntax), to keep one file under GitHub's 100 MB limit; the miRNA benchmark is split that way.
#
# Reads <results_root>/<method>/<case>/sample_<n>/rep<k>/summary.rds for every method directory under the root
# and writes <out_dir>/<dataset>_features.csv.gz with one row per method, case, sample size, replicate and
# feature, UNFILTERED: every feature a method's summary carries, with its `tested` flag where the method declined
# some, and no restriction to a common set. Columns are the union of the summaries' columns (NBSR: post_mean,
# post_sd, p_pos and eleven posterior quantiles of the shift-adjusted effect, centring_fraction; competitors:
# estimate, se, p_value, q_value, called, passed_sens, target_scale), plus `method`, the truth columns of fc.csv
# (alpha_null, alpha_alt, true_effect, is_de, ...) and, when a wald_comparison directory is given, NBSR's
# delta-method statistic (dm_est, dm_se, p_dm) joined by feature. Also writes <dataset>_manifest.csv (rows and
# units per method) and the wall time per unit from fit_complete.json in <dataset>_units.csv.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
method_re <- sub("^--methods=", "", grep("^--methods=", args, value = TRUE)); args <- grep("^--methods=", args, value = TRUE, invert = TRUE)
if (!length(method_re)) method_re <- "."
if (length(args) < 3) stop("Usage: export_paper_data.R <out_dir> <dataset> <results_root> [<wald_comparison_dir>]", call. = FALSE)
out_dir <- args[1]; dataset <- args[2]; root <- args[3]; wald_dir <- if (length(args) >= 4) args[4] else NA
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
methods <- basename(Filter(function(d) length(Sys.glob(file.path(d, "*", "sample_*", "rep*", "summary.rds"))) > 0,
                           list.dirs(root, recursive = FALSE)))
methods <- methods[grepl(method_re, methods, perl = TRUE)]
cat(dataset, ":", length(methods), "methods under", root, "\n")

units <- list(); D <- rbindlist(lapply(methods, function(m) {
  fs <- Sys.glob(file.path(root, m, "*", "sample_*", "rep*", "summary.rds"))
  d <- rbindlist(lapply(fs, function(f) {
    s <- as.data.table(readRDS(f))
    fc <- file.path(dirname(f), "fit_complete.json")
    if (file.exists(fc)) { j <- read_json(fc, simplifyVector = TRUE)
      units[[length(units) + 1]] <<- data.table(method = m, simul_case = s$simul_case[1], sample_size = s$sample_size[1],
        replicate_id = s$replicate_id[1], wall_seconds = j$total_seconds %||% NA_real_, chains = j$chains %||% NA_integer_,
        completed_at = j$completed_at %||% NA_character_) }
    s
  }), fill = TRUE)
  d[, method := m]
  if (!is.na(wald_dir)) {
    wf <- Sys.glob(file.path(wald_dir, paste0(m, "_by_feature.csv")))
    if (length(wf) == 1) {
      W <- fread(wf)[, .(simul_case, sample_size, replicate_id, miRNA, dm_est, dm_se, p_dm)]
      d <- merge(d, W, by = c("simul_case", "sample_size", "replicate_id", "miRNA"), all.x = TRUE)
      cat("  ", m, ": delta-method statistic joined from", basename(wf), "\n") } }
  d
}), fill = TRUE)
setnames(D, "miRNA", "feature")
if ("method_id" %in% names(D)) D[, method_id := NULL]
first <- intersect(c("method", "simul_case", "sample_size", "replicate_id", "feature", "stratum", "alpha_null", "alpha_alt",
                     "true_effect", "true_effect_abs", "is_de", "tested", "called", "passed_sens", "estimate", "se", "p_value", "q_value",
                     "post_mean", "post_sd", "p_pos", "dm_est", "dm_se", "p_dm"), names(D))
setcolorder(D, c(first, setdiff(names(D), first)))
setorder(D, method, simul_case, sample_size, replicate_id)
fwrite(D, file.path(out_dir, paste0(dataset, "_features.csv.gz")))
fwrite(rbindlist(units), file.path(out_dir, paste0(dataset, "_units.csv")))
man <- D[, .(units = uniqueN(paste(simul_case, sample_size, replicate_id)), rows = .N,
             cases = paste(sort(unique(simul_case)), collapse = " "), sample_sizes = paste(sort(unique(sample_size)), collapse = " ")), by = method]
fwrite(man, file.path(out_dir, paste0(dataset, "_manifest.csv"))); print(man)
cat("columns:", paste(names(D), collapse = ", "), "\n")
