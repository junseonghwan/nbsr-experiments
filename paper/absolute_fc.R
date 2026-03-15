rm(list=ls())

set.seed(1)

library(data.table)
library(DESeq2)
library(ggplot2)
library(MCMCpack)
library(seqinr)
library(stringr)

source("paper/dirichlet_estimation_functions.R")
source("paper/functions.R")

#output_path <- commandArgs(TRUE)[1]
sample_counts <- c(3, 5, 10, 20)
swap_no <- 2
output_path <- paste0("data/test/swap", swap_no)
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = T)
}
sample_count <- sample_counts[swap_no+1]
zero_miRNA_count <- 200
rep_count <- 20
swap_count <- 40
seq_depth <- 10^(6:8)

miRNA_polya_est <- fread(file = "data/microRNAome_dir_params.csv")

# Select features from mirGeneDB.
mirGeneDB <- read.fasta(file = "data/hsa.fas")
miRNA_names <- names(mirGeneDB)

lowercase_string <- tolower(miRNA_names)
corrected_case_string <- sub("mir", "miR", lowercase_string)
mirGeneDB_final_string <- gsub("_", "-", corrected_case_string)

pattern <- "hsa-(let|miR)-[1-9][0-9]*"
mirGeneDB_prefixes <- unlist(regmatches(mirGeneDB_final_string, gregexpr(pattern, mirGeneDB_final_string)))

# Find the prefix for each of the miRNAs.
miRNA_polya_est[,miRNA_prefix:=unlist(regmatches(miRNA, gregexpr(pattern, miRNA)))]
miRNA_polya_est <- miRNA_polya_est[miRNA_prefix %in% mirGeneDB_prefixes]

miRNA_polya_est[,sum(alpha_hat),by=.(cell_type)]
celltypes <- unique(miRNA_polya_est$cell_type)

unique(miRNA_polya_est$miRNA)

# Let's start simple.
ct <- celltypes[1]
alpha_hat <- miRNA_polya_est[cell_type == ct & alpha_hat > 0]
setorder(alpha_hat, -alpha_hat)
miRNA_count <- dim(alpha_hat)[1]

alpha_null <- copy(alpha_hat)
alpha_alt <- copy(alpha_null)

# Sample total abundance.
S0 <- sum(alpha_null$alpha_hat)
S1 <- sum(alpha_alt$alpha_hat)
true_bias <- -log2(S0/S1)

fc <- rep(1, miRNA_count)
K0 <- 20
# Sample log2fc
log2_fc <- rnorm(K0, mean = 0, sd = 3) 
fc[1:K0] <- 2^log2_fc

alpha_alt[,alpha_hat := fc * alpha_hat]

alpha_null[,sum(alpha_hat)]
alpha_alt[,sum(alpha_hat)]

alpha_null[,alpha_bar := normalize(alpha_hat)]
alpha_alt[,alpha_bar := normalize(alpha_hat)]

alpha_alt[,alpha := alpha_hat]
alpha_null[,alpha := alpha_hat]

# Generate counts.
ret <- generate_data(alpha_null, alpha_alt, sample_count, sample_count, reads_min_max = c(seq_depth[1], tail(seq_depth,1)), max_percent_miRNA = 0.5)

write.table(ret$X, file = paste0(output_path, "/X.csv"), quote = F, col.names = T, row.names = F, sep=",")
write.table(ret$Y, file = paste0(output_path, "/Y.csv"), quote = F, col.names = T, row.names = T, sep=",")
fwrite(ret$log_fc, file = paste0(output_path, "/fc.csv"))


plot(ret$log_fc$absolute_fc[1:K0], fc[1:K0])

se <- SummarizedExperiment(assays = list(counts = ret$Y), colData = ret$X)
dds2 <- DESeqDataSet(se, ~ trt)
dds2 <- DESeq(dds2, fitType = "local")
res2 <- results(dds2, contrast = c("trt", "alt", "null"))
deseq2_mu <- assays(dds2)[["mu"]]
deseq2_props <- apply(deseq2_mu, 2, normalize)

R <- colSums(ret$Y)
nbsr_pi <- read.csv(paste0(output_path, "/nbsr_pi.csv"), header = F)
nbsr_phi <- as.matrix(read.csv(paste0(output_path, "/nbsr_dispersion.csv"), header = F))
contrast_path <- paste0(output_path, "/alt_null")
nbsr_log2_fc <- log2(nbsr_pi[,sample_count+1]) - log2(nbsr_pi[,1])
nbsr_logRR <- as.matrix(read.csv(paste0(contrast_path, "/nbsr_logRR.csv"), header = F))
nbsr_logRR_sd <- as.matrix(read.csv(paste0(contrast_path, "/nbsr_logRR_sd.csv"), header = F))
nbsr_stats <- nbsr_logRR[1,] / nbsr_logRR_sd[1,]
nbsr_pvals <- 2*pnorm(abs(nbsr_stats), lower.tail = FALSE)
nbsr_adj_pvals <- p.adjust(nbsr_pvals, method = "BH")
nbsr_res <- data.table(features=rownames(ret$Y), log2FoldChange = nbsr_log2_fc, padj=nbsr_adj_pvals)

head(ret$log_fc)
which(res2$padj < 0.01)
which(nbsr_res$padj < 0.01)

# Biased by some constant: b.
# Compute true shift.
gt_bias <- -log2(sum(fc * alpha_null$alpha_bar))
print(abs(gt_bias - ret$log_fc$log2_fc[K0+1]))
plot(ret$log_fc$log2_fc, nbsr_res$log2FoldChange)
abline(a=0, b=1)
-log2(sum(fc * nbsr_pi[,1]))
-log2(sum(fc * deseq2_props[,1]))

not_sig_idx <- nbsr_res$padj < 0.01
# Find mode.
y <- nbsr_res$log2FoldChange[not_sig_idx]
median(y)
den <- density(y)
mode_idx <- which.max(den$y)
bias_est <- den$x[mode_idx]
print(den$x[mode_idx])
print(gt_bias)

temp <- data.table(gt=fc, deseq2=2^res2$log2FoldChange, nbsr=2^nbsr_log2_fc)
head(temp, K0)
temp[1:K0,sqrt(mean((gt - deseq2)^2))]
temp[1:K0,sqrt(mean((gt - nbsr)^2))]

debiased_fc_dt <- data.table(gt=fc, deseq2=2^res2$log2FoldChange, nbsr=2^(nbsr_log2_fc-bias_est))
head(debiased_fc_dt, K0)

# compute cumulative sum of square of errors.
deseq2_rmse <- sqrt(cumsum(debiased_fc_dt[,(deseq2 - gt)^2])/(1:miRNA_count))
nbsr_rmse <- sqrt(cumsum(debiased_fc_dt[,(nbsr - gt)^2])/(1:miRNA_count))

# Check correctness of RMSE computation.
temp2[1:K0,sqrt(mean((gt - deseq2)^2))]
temp2[1:K0,sqrt(mean((gt - nbsr)^2))]
deseq2_rmse[K0]
nbsr_rmse[K0]

plot(nbsr_rmse, deseq2_rmse) # Something weird  happens for DESeq2.
plot(1:miRNA_count, nbsr_rmse) # Errors are very controlled for NBSR.
