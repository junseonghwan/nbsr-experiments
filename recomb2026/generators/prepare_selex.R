# SELEX (McMurrough et al. 2014, PNAS; the ALDEx2 package's `selex`): 1600 I-LtrI variants, 7 non-selected and
# 7 selected samples. Writes the pipeline layout under <out>:
#   SELEX/sample_7/rep1         the full data
#   SELEX/sample_5/rep1..rep20  random 5 + 5 subsets without replacement (as Nixon, Gloor & Silverman 2025)
# fc.csv truth: the 27 variants McMurrough et al. characterised individually (`nixon27`), the reference set of
# Fernandes et al. (2014) and Nixon et al. (2025). Their individual measurements ride along per variant -- survival
# on selective medium in the I-LtrI and I-OnuI backbones, and minutes to mid-log growth -- without being used to
# split the set. The selection is bacteriostatic, so a variant that cannot clear the toxin stops dividing and stays
# at its input abundance: relative to that majority no variant can be down.
# Usage: Rscript data_processor/prepare_selex.R <out> [--subsets 20] [--seed 1]
args <- commandArgs(TRUE); out <- args[1]
opt <- function(f, d) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
n_sub <- as.integer(opt("--subsets", "20")); set.seed(as.integer(opt("--seed", "1")))
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
Y <- as.matrix(read.csv("data/selex/Y.csv", row.names = 1, check.names = FALSE))
X <- read.csv("data/selex/X.csv", row.names = 1, check.names = FALSE)
truth <- read.delim("data/selex/truth.txt", row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(Y), rownames(X)), all(rownames(truth) %in% rownames(Y)))
characterised <- rownames(Y) %in% rownames(truth)
fc <- data.frame(miRNA = rownames(Y), alpha_null = 1, alpha_alt = ifelse(characterised, exp(1), 1),
                 nixon27 = characterised, ltr_survival = truth[rownames(Y), "LTR"],
                 onu_survival = truth[rownames(Y), "ONU"], minutes_to_midlog = truth[rownames(Y), "T"])
fc$stratum <- "all"
ns <- which(X$trt == "NS"); s <- which(X$trt == "S")
write_unit <- function(n, rep, cols) {
  d <- file.path(out, "SELEX", paste0("sample_", n), paste0("rep", rep)); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.csv(Y[, cols], file.path(d, "Y.csv"))
  write.csv(data.frame(trt = ifelse(X$trt[cols] == "S", "alt", "null"), row.names = colnames(Y)[cols]), file.path(d, "X.csv"))
  write.csv(fc, file.path(d, "fc.csv"), row.names = FALSE)
}
write_unit(7, 1, c(ns, s))
for (r in seq_len(n_sub)) write_unit(5, r, c(sort(sample(ns, 5)), sort(sample(s, 5))))
cat("SELEX written to", out, ":", sum(characterised), "characterised variants,", sum(!characterised), "others; 1 full unit and", n_sub, "subsets of 5 + 5\n")
