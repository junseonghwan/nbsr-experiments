# Label-permutation null from real microbiome data: the two arms of a real comparison are pooled and the labels
# reassigned at random, so every genus is null with respect to the label while the real between-subject and
# between-group variation stays in the counts. A method's calls on these replicates are false positives.
#
# Usage: Rscript data_processor/prepare_permutation.R <out> [--sizes 20,40,80] [--reps 50] [--seed 1]
#
# Datasets (data/microbiome/processed/<template>): the case name is the pipeline's "simulation case", so one data
# root holds both and a single node fits them together.
#   BAXTER    crc_baxter_genus,   H (172) + CRC (120)
#   SCHUBERT  cdi_schubert_genus, H (154) + CDI (93)
# Each replicate draws 2n samples from the pool without replacement and assigns n to each arm, so the mix of
# healthy and disease samples in an arm varies at random; that heterogeneity is what the test exposes. The seed is
# base_seed * 1e7 + case_id * 1e5 + n * 1e2 + rep, fixed per case, so adding a case or a size leaves the rest
# unchanged and every method sees identical replicates.
#
# Output (the layout the pipeline reads): <out>/<CASE>/sample_<n>/rep<k>/{Y.csv, X.csv, fc.csv, settings.json}
# plus <out>/replicate_summary.csv. fc.csv gives every genus alpha_null = alpha_alt = 1, i.e. true_effect = 0 and
# is_de = FALSE, and carries the stratum by mean proportion for reporting. settings.json records which original
# group each arm's samples came from.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
args <- commandArgs(TRUE)
opt <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
out_root <- args[1]
if (is.na(out_root)) stop("Usage: prepare_permutation.R <out> [--sizes 20,40,80] [--reps 50] [--seed 1]", call. = FALSE)
sizes     <- as.integer(strsplit(opt("--sizes", "20,40,80"), ",")[[1]])
n_reps    <- as.integer(opt("--reps", "50"))
base_seed <- as.integer(opt("--seed", "1"))

CASES <- list(
  BAXTER   = list(seed_id = 1L, template = "crc_baxter_genus",   groups = c("H", "CRC")),
  SCHUBERT = list(seed_id = 2L, template = "cdi_schubert_genus", groups = c("H", "CDI")))
STRATA_CUTS <- c(HIGH = 0.01, MED = 0.001)   # mean proportion, as in the simulator
MIN_SAMPLES <- 1                             # drop a genus only if it has no read in any sample of the replicate

summary_rows <- list()
for (case in names(CASES)) {
  cfg <- CASES[[case]]
  d <- file.path("data/microbiome/processed", cfg$template)
  Y <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE))
  X <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
  stopifnot(identical(colnames(Y), rownames(X)), all(cfg$groups %in% X$trt))
  pool <- which(X$trt %in% cfg$groups)
  cat(sprintf("%s (%s): pool %d samples (%s), %d genera\n", case, cfg$template, length(pool),
              paste(sprintf("%s=%d", cfg$groups, sapply(cfg$groups, function(g) sum(X$trt[pool] == g))), collapse = " "),
              nrow(Y)))
  for (n in sizes) {
    if (2 * n > length(pool)) stop("case ", case, ": 2n = ", 2 * n, " exceeds the pool of ", length(pool), call. = FALSE)
    for (rep in seq_len(n_reps)) {
      set.seed(base_seed * 1e7 + cfg$seed_id * 1e5 + n * 1e2 + rep)
      picked <- sample(pool, 2 * n)                       # without replacement within a replicate
      arm <- rep(c("null", "alt"), each = n)   # the pipeline's control and treatment levels
      Yc <- Y[, picked, drop = FALSE]
      retained <- rowSums(Yc > 0) >= MIN_SAMPLES
      Yc <- Yc[retained, , drop = FALSE]
      prop <- rowMeans(sweep(Yc, 2, colSums(Yc), "/"))
      fc <- data.table(miRNA = rownames(Yc), alpha_null = 1, alpha_alt = 1,
                       stratum = fifelse(prop > STRATA_CUTS["HIGH"], "HIGH",
                                         fifelse(prop > STRATA_CUTS["MED"], "MED", "LOW")),
                       mean_prop = prop)
      dir_out <- file.path(out_root, case, sprintf("sample_%d", n), sprintf("rep%d", rep))
      dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
      write.csv(Yc, file.path(dir_out, "Y.csv"))
      write.csv(data.frame(trt = arm, row.names = colnames(Yc)), file.path(dir_out, "X.csv"))
      fwrite(fc, file.path(dir_out, "fc.csv"))
      orig <- X$trt[picked]
      write_json(list(case = case, template = cfg$template, pooled_groups = cfg$groups, sample_size = n,
                      replicate = rep, seed = base_seed * 1e7 + cfg$seed_id * 1e5 + n * 1e2 + rep,
                      arm_composition = list(null = as.list(table(orig[arm == "null"])),
                                             alt  = as.list(table(orig[arm == "alt"]))),
                      n_genera_template = nrow(Y), n_genera_retained = nrow(Yc)),
                 file.path(dir_out, "settings.json"), auto_unbox = TRUE, pretty = TRUE)
      summary_rows[[length(summary_rows) + 1]] <- data.table(
        case = case, sample_size = n, replicate = rep, n_retained = nrow(Yc),
        zero_frac = mean(Yc == 0), lib_median = median(colSums(Yc)),
        disease_in_null = sum(orig[arm == "null"] == cfg$groups[2]),
        disease_in_alt  = sum(orig[arm == "alt"]  == cfg$groups[2]))
    }
  }
}
summ <- rbindlist(summary_rows)
fwrite(summ, file.path(out_root, "replicate_summary.csv"))
cat("\nWritten to", out_root, ":", nrow(summ), "replicates\n")
print(summ[, .(reps = .N, genera = round(mean(n_retained)), zero_frac = round(mean(zero_frac), 2),
               disease_share_null = round(mean(disease_in_null / sample_size), 2),
               imbalance_q90 = round(quantile(abs(disease_in_null - disease_in_alt) / sample_size, 0.9), 2)),
         by = .(case, sample_size)][order(case, sample_size)])
