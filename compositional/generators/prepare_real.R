# Real two-group analyses with the original labels, from the processed MicrobiomeHD studies the label-permutation
# null (prepare_permutation.R) draws on: Baxter et al. (2016) stool, healthy against carcinoma and against adenoma;
# Schubert et al. (2014) stool, healthy against C. difficile infection.
#
# Usage: Rscript data_processor/prepare_real.R <out> [CASE ...]
#
# Output in the pipeline's layout, one replicate per case: <out>/<CASE>/sample_<n_control>/rep1/{Y.csv, X.csv,
# fc.csv, settings.json}. The size label is the control arm, so the launcher is called per dataset with that size.
# The control level is "null" (healthy) and the treatment level "alt" (disease), so a positive effect is enrichment
# in disease. fc.csv carries no truth (alpha_null = alpha_alt = 1) and the stratum by mean proportion.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
args <- commandArgs(TRUE); out_root <- args[1]
if (is.na(out_root)) stop("Usage: prepare_real.R <out> [CASE ...]", call. = FALSE)
CASES <- list(BAXTER_CRC     = list(template = "crc_baxter_genus",   groups = c("H", "CRC")),
              BAXTER_ADENOMA = list(template = "crc_baxter_genus",   groups = c("H", "nonCRC")),
              SCHUBERT_CDI   = list(template = "cdi_schubert_genus", groups = c("H", "CDI")))
wanted <- if (length(args) > 1) args[-1] else names(CASES)
STRATA_CUTS <- c(HIGH = 0.01, MED = 0.001)
for (case in wanted) {
  cfg <- CASES[[case]]; g <- cfg$groups; d <- file.path("data/microbiome/processed", cfg$template)
  Y <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE))
  X <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
  stopifnot(identical(colnames(Y), rownames(X)), all(g %in% X$trt))
  keep <- which(X$trt %in% g)
  Yc <- Y[, keep, drop = FALSE]; Yc <- Yc[rowSums(Yc > 0) >= 1, , drop = FALSE]
  arm <- ifelse(X$trt[keep] == g[1], "null", "alt")
  prop <- rowMeans(sweep(Yc, 2, colSums(Yc), "/"))
  fc <- data.table(miRNA = rownames(Yc), alpha_null = 1, alpha_alt = 1,
                   stratum = fifelse(prop > STRATA_CUTS["HIGH"], "HIGH", fifelse(prop > STRATA_CUTS["MED"], "MED", "LOW")),
                   mean_prop = prop)
  n_ctrl <- sum(arm == "null")
  dir_out <- file.path(out_root, case, sprintf("sample_%d", n_ctrl), "rep1")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  write.csv(Yc, file.path(dir_out, "Y.csv"))
  write.csv(data.frame(trt = arm, row.names = colnames(Yc)), file.path(dir_out, "X.csv"))
  fwrite(fc, file.path(dir_out, "fc.csv"))
  write_json(list(case = case, template = cfg$template, groups = list(null = g[1], alt = g[2]),
                  n_null = n_ctrl, n_alt = sum(arm == "alt"), n_genera = nrow(Yc), zero_frac = mean(Yc == 0),
                  lib_median = median(colSums(Yc))),
             file.path(dir_out, "settings.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("%s: %s=%d vs %s=%d, %d genera, zero fraction %.2f, median depth %d -> %s\n", case, g[1], n_ctrl, g[2],
              sum(arm == "alt"), nrow(Yc), mean(Yc == 0), as.integer(median(colSums(Yc))), dir_out))
}
