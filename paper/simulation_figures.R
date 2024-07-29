# Generate figure for each simulation scenario.
rm(list=ls())
library(cowplot)
library(data.table)
library(dplyr)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(stringr)
#library(rlist)
source("paper/functions.R")

swap_no <- 0
USE_EB <- ifelse(swap_no == 2, TRUE, FALSE) # for swap_no=2, we will process swap2/ swap2_EB/ together.
data_path <- paste0("data/validation/swap", swap_no)
results <- fread(paste0(data_path, "/results.csv"))
performance <- fread(paste0(data_path, "/performance.csv"))

# NBSR with fixed dispersion.
if (USE_EB) {
  data_path2 <- paste0(data_path, "_EB/")
  results2 <- fread(paste0(data_path2, "results.csv"))
  performance2 <- fread(paste0(data_path2, "performance.csv"))
  
  results2[method == "NBSR",method := "NBSR EB"]
  performance2[method == "NBSR",method := "NBSR EB"]
  results <- rbind(results, results2[method == "NBSR EB"])
  performance <- rbind(performance, performance2[method == "NBSR EB"])
}

results$celltype2 <- gsub(pattern = "_", replacement = " ", x = stringr::str_to_title(results$celltype))
results$celltype2 <- gsub(pattern = "cd", replacement = "CD", results$celltype2)
results$celltype2 <- gsub(pattern = " CD19| CD56", replacement = "", results$celltype2)

performance$celltype2 <- gsub(pattern = "_", replacement = " ", x = stringr::str_to_title(performance$celltype))
performance$celltype2 <- gsub(pattern = "cd", replacement = "CD", performance$celltype2)
performance$celltype2 <- gsub(pattern = " CD19| CD56", replacement = "", performance$celltype2)

method_levels <- c("DESeq2", "EdgeR", "NBSR EB", "NBSR")
results[,method := factor(method, levels=method_levels)]
performance[,method := factor(method, levels=method_levels)]

#figure_path <- "../NBSR/paper/figures/"
figure_path <- paste0("paper/figures/swap", swap_no, "/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = TRUE)
}

# Generate RMSE and coverage plots.
cols <- RColorBrewer::brewer.pal(8, "Set1")
cols2 <- RColorBrewer::brewer.pal(8, "Set2")
pl <- ggplot(results, aes(method, rmse, fill=method)) + geom_boxplot()
pl <- pl + theme_cowplot()
pl <- pl + ylab("RMSE") + xlab("") + theme(legend.title = element_blank())
pl <- pl + facet_grid(~ celltype2)
pl <- decorate_figure(pl, xlab_text = "", ylab_text = "RMSE")
pl <- pl + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pl <- pl + scale_fill_manual(values=cols[1:4])
ggsave(filename = paste0(figure_path, "/rmse.pdf"), pl, width=14, height = 4)

pl <- ggplot(results[method != "EdgeR"], aes(method, coverage, fill=method)) + geom_boxplot()
pl <- pl + theme_cowplot()
pl <- pl + ylab("Coverage") + xlab("") + theme(legend.title = element_blank())
#pl <- pl + scale_fill_manual(values=cols[1:2])
pl <- pl + geom_hline(yintercept = 0.95, linetype=2)
pl <- pl + facet_grid(~ celltype2)
pl <- decorate_figure(pl, xlab_text = "", ylab_text = "Coverage")
pl <- pl + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pl <- pl + scale_fill_manual(values=cols[c(1:2, 4)])
ggsave(filename = paste0(figure_path, "/coverage.pdf"), pl, width = 14)

# Generate precision and false positive rate vs threshold.
precision_dt <- performance[,.(avg=mean(prec), std=sd(prec)),by=.(threshold, method)]
pl <- ggplot(precision_dt, aes(factor(threshold), avg)) + geom_point() + geom_line(aes(group=method))
pl <- pl + theme_cowplot()
pl <- pl + geom_errorbar(aes(ymin = avg-2*std, ymax = avg+2*std), width=0.1)
#pl <- pl + scale_color_manual(values=cols[1:3])
pl <- pl + ylab("Precision") + xlab("Significance threshold") + theme(legend.title = element_blank())
pl <- pl + facet_grid(~method)
ggsave(filename = paste0(figure_path, "/precision.pdf"), pl)

fpr <- performance[,.(avg=mean(fpr), std=sd(fpr)),by=.(threshold, method)]
fnr <- performance[,.(avg=mean(fnr), std=sd(fnr)),by=.(threshold, method)]
fpr$type <- "FPR"
fnr$type <- "FNR"

pl <- ggplot(fpr, aes(factor(threshold), avg)) + geom_point() + geom_line(aes(group=method))
pl <- pl + theme_cowplot()
pl <- pl + geom_errorbar(aes(ymin = avg-2*std, ymax = avg+2*std), width=0.1)
#pl <- pl + scale_color_manual(values=cols[1:3])
pl <- pl + ylab("FPR") + xlab("Significance threshold") + theme(legend.title = element_blank())
pl <- pl + facet_grid(~ method)
ggsave(filename = paste0(figure_path, "/fpr.pdf"), pl)

pl <- ggplot(fnr, aes(factor(threshold), avg)) + geom_point() + geom_line(aes(group=method))
pl <- pl + theme_cowplot(12)
pl <- pl + geom_errorbar(aes(ymin = avg-2*std, ymax = avg+2*std), width=0.1)
#pl <- pl + scale_color_manual(values=cols[1:3])
pl <- pl + ylab("FNR") + xlab("Significance threshold") + theme(legend.title = element_blank())
pl <- pl + facet_grid(~ method)
ggsave(filename = paste0(figure_path, "/fnr.pdf"), pl)

# Plot FPR vs FNR as threshold is varied.
error_rates <- rbind(fpr, fnr)
pl <- ggplot(error_rates, aes(factor(threshold), avg, fill=type)) + geom_bar(stat = "identity", position = "stack")
pl <- pl + facet_grid(~ method)
pl <- pl + theme_cowplot()
pl <- pl + scale_fill_manual(values=cols2[1:2])
pl <- pl + xlab("Significance threshold") + ylab("Error rate") + theme(legend.title = element_blank())
ggsave(filename = paste0(figure_path, "/fpr_fnr.pdf"), pl)

# Dispersion plot.
if (swap_no >= 2) {
  
  # Select cell type and rep.
  # Compute true biological CV.
  # Compare estimated dispersion.
  celltypes <- unique(results$celltype)
  celltype <- celltypes[1]
  rep_no <- 1
  rep_path <- paste0(data_path, "/", celltype, "/rep", rep_no)
  Y <- read.csv(paste0(rep_path, "/Y.csv"))
  X <- read.csv(paste0(rep_path, "/X.csv"))
  log_fc <- fread(paste0(rep_path, "/fc.csv"))
  dat_null <- data.table(alpha=log_fc$alpha_null, alpha_bar=log_fc$alpha_null_bar)
  dat_alt <- data.table(alpha=log_fc$alpha_alt, alpha_bar=log_fc$alpha_alt_bar)
  compute_cv(dat_null)
  compute_cv(dat_alt)
  log_fc <- log_fc[miRNA %in% rownames(Y),]
  feature_count <- dim(Y)[1]
  sample_count <- dim(Y)[2]
  
  feature_sparsity <- rowMeans(Y == 0)
  
  Y <- as.matrix(Y)
  R <- colSums(Y)
  
  # DESeq2
  se <- SummarizedExperiment(assays = list(counts = Y), colData = X)
  dds2 <- DESeqDataSet(se, ~ trt)
  dds2 <- DESeq(dds2)
  res2 <- results(dds2, contrast = c("trt", "alt", "null"))
  deseq2_phi <- rowData(dds2)$dispersion
  sum(is.na(deseq2_phi))
  deseq2_phi[is.na(deseq2_phi)] <- 0
  mean_expr <- assays(dds2)[["mu"]]
  rdat <- rowData(dds2)
  
  # EdgeR
  d <- edgeR::DGEList(counts=Y)
  d <- calcNormFactors(d)
  design_mat <- model.matrix(~ trt, X)
  d <- edgeR::estimateDisp(d, design = design_mat)
  fit <- glmFit(d, design_mat)
  lrt <- glmLRT(fit, contrast=c(0, -1))
  edgeR_results <- topTags(lrt, n=Inf, sort.by="none")
  edgeR_phi <- fit$dispersion
  edgeR_mu <- fit$fitted.values
  
  # NBSR -- load the results.
  nbsr_pi <- read.csv(paste0(rep_path, "/nbsr_pi.csv"), header = F)
  nbsr_phi <- as.matrix(read.csv(paste0(rep_path, "/nbsr_dispersion.csv"), header = F))
  nbsr_mu <- sweep(nbsr_pi, 2, R, "*")
  
  rownames(nbsr_pi) <- rownames(nbsr_mu) <- rownames(Y)
  colnames(nbsr_pi) <- colnames(nbsr_mu) <- colnames(Y)
  
  nbsr_pi <- as.matrix(nbsr_pi)
  nbsr_mu <- as.matrix(nbsr_mu)
  
  ## Plot dispersion estimates for perturbed features.
  set.seed(1)
  gt_phi <- cbind(dat_null$phi, dat_alt$phi)
  
  null_idx <- X$trt == "null"
  alt_idx <- !null_idx
  
  perturb_idxs <- log_fc$log_fc != 0
  
  nbsr_phi_null <- nbsr_phi[perturb_idxs,null_idx]
  nbsr_phi_alt <- nbsr_phi[perturb_idxs,alt_idx]
  
  # NBSR EB
  if (USE_EB) {
    rep_path2 <- paste0(data_path2, "/", celltype, "/rep", rep_no)
    nbsr_eb_pi <- read.csv(paste0(rep_path2, "/nbsr_pi.csv"), header = F)
    nbsr_eb_phi <- as.matrix(read.csv(paste0(rep_path2, "/nbsr_dispersion.csv"), header = F))[,1]
    nbsr_eb_mu <- sweep(nbsr_eb_pi, 2, R, "*")
    
    rownames(nbsr_eb_pi) <- rownames(nbsr_eb_mu) <- rownames(Y)
    colnames(nbsr_eb_pi) <- colnames(nbsr_eb_mu) <- colnames(Y)
    
    nbsr_eb_pi <- as.matrix(nbsr_eb_pi)
    nbsr_eb_mu <- as.matrix(nbsr_eb_mu)
    
    null_phi_dt <- data.table(feature = rownames(Y)[perturb_idxs], 
                              gt = gt_phi[perturb_idxs,1],
                              edgeR_phi = edgeR_phi[perturb_idxs],
                              deseq2_phi = deseq2_phi[perturb_idxs],
                              nbsr_eb_phi = nbsr_eb_phi[perturb_idxs],
                              nbsr_phi = nbsr_phi_null)
    null_phi_dt <- melt(null_phi_dt, id.vars = c("feature", "gt", "deseq2_phi", "edgeR_phi", "nbsr_eb_phi"), variable.name = "sample", value.name = "nbsr_phi")
    alt_phi_dt <- data.table(feature = rownames(Y)[perturb_idxs], 
                             gt = gt_phi[perturb_idxs,2],
                             deseq2_phi = deseq2_phi[perturb_idxs],
                             edgeR_phi = edgeR_phi[perturb_idxs],
                             nbsr_eb_phi = nbsr_eb_phi[perturb_idxs],
                             nbsr_phi = nbsr_phi_alt)
    alt_phi_dt <- melt(alt_phi_dt, id.vars = c("feature", "gt", "deseq2_phi", "edgeR_phi", "nbsr_eb_phi"), variable.name = "sample", value.name = "nbsr_phi")
  } else {
    null_phi_dt <- data.table(feature = rownames(Y)[perturb_idxs], 
                              gt = gt_phi[perturb_idxs,1],
                              edgeR_phi = edgeR_phi[perturb_idxs],
                              deseq2_phi = deseq2_phi[perturb_idxs],
                              nbsr_phi = nbsr_phi_null)
    null_phi_dt <- melt(null_phi_dt, id.vars = c("feature", "gt", "deseq2_phi", "edgeR_phi"), variable.name = "sample", value.name = "nbsr_phi")
    alt_phi_dt <- data.table(feature = rownames(Y)[perturb_idxs], 
                             gt = gt_phi[perturb_idxs,2],
                             deseq2_phi = deseq2_phi[perturb_idxs],
                             edgeR_phi = edgeR_phi[perturb_idxs],
                             nbsr_phi = nbsr_phi_alt)
    alt_phi_dt <- melt(alt_phi_dt, id.vars = c("feature", "gt", "deseq2_phi", "edgeR_phi"), variable.name = "sample", value.name = "nbsr_phi")
  }
  
  null_phi_dt$trt <- "null"
  alt_phi_dt$trt <- "alt"
  phi_dt <- rbind(null_phi_dt, alt_phi_dt)
  
  feature_idx <- sample(which(perturb_idxs), 5)
  fnames <- rownames(Y)[feature_idx]
  sub_dt <- phi_dt[feature %in% fnames]
  sub_dt[,feature:=gsub("hsa-", "", feature)]
  
  pl <- ggplot(sub_dt, aes(trt, nbsr_phi)) + theme_cowplot() +
    geom_boxplot(aes(fill="NBSR")) +
    geom_point(aes(trt, gt, shape="True dispersion"), size=5) +  
    geom_segment(aes(x = -Inf, xend = Inf, y = deseq2_phi, yend = deseq2_phi, color = 'DESeq2'), 
                 linetype = "dashed", linewidth = 1) +  
    geom_segment(aes(x = -Inf, xend = Inf, y = edgeR_phi, yend = edgeR_phi, color = 'EdgeR'), 
                 linetype = "dotted", linewidth = 1) +
    xlab("Experimental condition") + 
    ylab("Biological CV") +
    facet_grid(~ feature) +
    scale_fill_manual(name="", values = c("NBSR" = cols[4])) + 
    scale_y_log10()
  if (USE_EB) {
    pl <- pl +   geom_segment(aes(x = -Inf, xend = Inf, y = nbsr_eb_phi, yend = nbsr_eb_phi, color = 'NBSR_EB'), 
                              linetype = "solid", linewidth = 1) +
      scale_color_manual(name = "", 
                         values = c("DESeq2" = cols[1], "EdgeR" = cols[2], "NBSR_EB" = cols[3]), 
                         labels = c("DESeq2" = "DESeq2", "EdgeR" = "EdgeR", "NBSR_EB"="NBSR EB")) +
      scale_shape_manual(name = "",
                         values = c("True dispersion" = 4),
                         labels = c("True dispersion" = "True dispersion"))
  } else {
    pl <- pl + scale_color_manual(name = "", 
                                  values = c("DESeq2" = cols[1], "EdgeR" = cols[2]), 
                                  labels = c("DESeq2" = "DESeq2", "EdgeR" = "EdgeR")) +
      scale_shape_manual(name = "",
                         values = c("True dispersion" = 4),
                         labels = c("True dispersion" = "True dispersion"))
  }
  pl <- decorate_figure(pl, ylab="Dispersion (Log10)", xlab="")
  ggsave(filename = paste0(figure_path, "/biological_cv.pdf"), pl, width=12)

}  

