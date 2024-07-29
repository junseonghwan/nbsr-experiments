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
swap_no <- 0
output_path <- paste0("data/validation/swap", swap_no)
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

# Swap parameter values across condition A vs B.
# Bin the miRNAs into groups: LOW, MED, HIGH.
# Swap LOW<->MED and LOW<->HIGH.
# Swapping within LOW bin will not yield any significant changes.
#miRNA_polya_est[,mean(0 <= alpha_hat & alpha_hat < 1),by=.(cell_type)]
#miRNA_polya_est[,mean(1 <= alpha_hat & alpha_hat < 10),by=.(cell_type)]
#miRNA_polya_est[,mean(10 <= alpha_hat & alpha_hat < Inf),by=.(cell_type)]

for (celltype in celltypes)
{
  alpha_hat <- miRNA_polya_est[cell_type == celltype & alpha_hat > 0]
  miRNA_count <- dim(alpha_hat)[1]
  
  # Set a subset of the miRNA expression matrix to 0 so that they are not expressed by any sample.
  for (rep_no in 1:rep_count)
  {
    alpha_null <- copy(alpha_hat)
    alpha_alt <- copy(alpha_null)
    
    low_idxs <- which(alpha_null$alpha_hat <= 1)
    high_idxs <- which(alpha_null$alpha_hat > 1)
    length(high_idxs)
    
    # Sample pairs of indices
    idx1 <-  sample(x = low_idxs, size = swap_count, replace = FALSE)
    idx2 <-  sample(x = high_idxs, size = swap_count, replace = FALSE)
    swap_idx_pairs <- cbind(idx1, idx2)
    
    alpha_alt$alpha_hat[idx1] <- alpha_null$alpha_hat[idx2]
    alpha_alt$alpha_hat[idx2] <- alpha_null$alpha_hat[idx1]
    
    alpha_null[,alpha_bar := normalize(alpha_hat)]
    alpha_alt[,alpha_bar := normalize(alpha_hat)]
    
    alpha_alt[,alpha := alpha_hat]
    alpha_null[,alpha := alpha_hat]
    
    low_idx <- which(alpha_null$alpha_hat < 1 & alpha_alt$alpha_hat < 1)
    zero_out_idx <- sample(low_idx, zero_miRNA_count, replace = F)
    alpha_alt$alpha[zero_out_idx] <- alpha_null$alpha[zero_out_idx] <- 1e-3
    
    ret <- generate_data(alpha_alt, alpha_null, sample_count, sample_count, reads_min_max = c(seq_depth[1], tail(seq_depth,1)), max_percent_miRNA = 0.5)
    
    # Output data to run NBSR.
    rep_path <- paste0(output_path, "/", celltype, "/rep", rep_no)
      dir.create(rep_path, recursive = TRUE)
    if (!dir.exists(rep_path)) {
    }

    row_idxs <- rowSums(ret$Y > 0) > 1
    sum(row_idxs)
    Y <- ret$Y[row_idxs,]
    fc <- ret$log_fc[row_idxs]
    write.table(ret$X, file = paste0(rep_path, "/X.csv"), quote = F, col.names = T, row.names = F, sep=",")
    write.table(Y, file = paste0(rep_path, "/Y.csv"), quote = F, col.names = T, row.names = T, sep=",")
    fwrite(fc, file = paste0(rep_path, "/fc.csv"))
  }
}
