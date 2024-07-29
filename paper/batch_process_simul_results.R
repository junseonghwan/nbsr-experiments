# getwd() points to nbsr-experiments/
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

miRNA_polya_est <- fread(file = "data/microRNAome_dir_params.csv")
celltypes <- unique(miRNA_polya_est$cell_type)
swap_no <- 0
data_path <- paste0("data/validation/swap", swap_no)
dispersion_model <- ifelse(swap_no < 2, "eb", "regression")
rep_count <- 20
thresholds <- c(0.01, 0.05, 0.1)

results_list <- list()
performance_list <- list()
for (celltype in celltypes)
{
  for (rep_no in 1:rep_count)
  {
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
    
    # NBSR
    nbsr_pi <- read.csv(paste0(rep_path, "/nbsr_pi.csv"), header = F)
    nbsr_phi <- as.matrix(read.csv(paste0(rep_path, "/nbsr_dispersion.csv"), header = F))
    nbsr_mu <- sweep(nbsr_pi, 2, R, "*")
    
    rownames(nbsr_pi) <- rownames(nbsr_mu) <- rownames(Y)
    colnames(nbsr_pi) <- colnames(nbsr_mu) <- colnames(Y)
    
    nbsr_pi <- as.matrix(nbsr_pi)
    nbsr_mu <- as.matrix(nbsr_mu)
    
    # Compute log2 fold change for NBSR.
    nbsr_log2_fc <- log2(nbsr_pi[,which(X$trt == "alt")[1]]) - log2(nbsr_pi[,which(X$trt == "null")[1]])

    deseq2_err <- log_fc$log2_fc - res2$log2FoldChange
    edgeR_err <- log_fc$log2_fc - edgeR_results$table$logFC
    nbsr_err <- log_fc$log2_fc - nbsr_log2_fc
    
    print(paste0("DESeq2 RMSE:",  sqrt(mean(deseq2_err^2))))
    print(paste0("EdgeR RMSE:",  sqrt(mean(edgeR_err^2))))
    print(paste0("NBSR RMSE:",  sqrt(mean(nbsr_err^2))))
    
    plt_dt <- data.table(deseq2=deseq2_err, edgeR=edgeR_err, nbsr=nbsr_err, log2fc=log_fc$log2_fc)
    plt_dt[,sign := ifelse(log2fc > 0, "+", ifelse(log2fc < 0, "-", "0"))]
    plt_dt$sign <- factor(plt_dt$sign, levels = c("-", "0", "+"))
    plt_dt <- melt(plt_dt, id.vars = c("sign", "log2fc"), variable.name = "method", value.name = "error")
    pl <- ggplot(plt_dt, aes(method, error, fill=sign)) + geom_boxplot()
    ggsave(filename = paste0(rep_path, "/error.png"), pl)
    
    deseq2_lb <- res2$log2FoldChange - 1.96*res2$lfcSE
    deseq2_ub <- res2$log2FoldChange + 1.96*res2$lfcSE
    deseq2_coverage <- (log_fc$log2_fc > deseq2_lb) & (log_fc$log2_fc < deseq2_ub)
    
    # EdgeR does not offer standard error estimates: https://support.bioconductor.org/p/61640/#61662
    
    nbsr_logRR <- as.matrix(read.csv(paste0(rep_path, "/nbsr_logRR.csv"), header = F))
    nbsr_logRR_sd <- as.matrix(read.csv(paste0(rep_path, "/nbsr_logRR_sd.csv"), header = F))
    
    nbsr_log_fc <- nbsr_logRR[1,]
    nbsr_log_fc_sd <- nbsr_logRR_sd[1,]
    nbsr_lb <- (nbsr_log_fc - 1.96*nbsr_log_fc_sd)
    nbsr_ub <- (nbsr_log_fc + 1.96*nbsr_log_fc_sd)
    nbsr_coverage <- (log_fc$log_fc > nbsr_lb) & (log_fc$log_fc < nbsr_ub)

    results_dt <- rbind(data.table(coverage=mean(deseq2_coverage, na.rm = T),
                                   rmse=sqrt(mean(deseq2_err^2, na.rm=T)),
                                   method="DESeq2"),
                        data.table(coverage=NA,
                                   rmse=sqrt(mean(edgeR_err^2, na.rm=T)),
                                   method="EdgeR"),
                        data.table(coverage=mean(nbsr_coverage, na.rm = T),
                                   rmse=sqrt(mean(nbsr_err^2, na.rm=T)),
                                   method="NBSR"))
    results_dt$celltype <- celltype
    results_dt$rep_no <- rep_no

    # Set the pvalues for features with NA to 1.
    res2$padj[is.na(res2$padj)] <- 1
    
    nbsr_stats <- nbsr_log_fc/nbsr_log_fc_sd
    nbsr_pvals <- 2*pnorm(abs(nbsr_stats), lower.tail = FALSE)
    nbsr_adj_pvals <- p.adjust(nbsr_pvals, method = "BH")
    
    gt_sig_idx <- which(log_fc$log_fc != 0)
    gt_not_sig_idx <- which(log_fc$log_fc == 0)

    deseq2_dt <- data.table(log_fc=res2$log2FoldChange, lb=deseq2_lb, ub=deseq2_ub, pval=res2$padj, gt_log_fc=log_fc$log2_fc, method="DESeq2")
    edgeR_dt <- data.table(log_fc=edgeR_results$table$logFC, lb=NA, ub=NA, pval=edgeR_results$table$FDR, gt_log_fc=log_fc$log2_fc, method="EdgeR")
    nbsr_dt <- data.table(log_fc=nbsr_log_fc, lb=nbsr_lb, ub=nbsr_ub, pval=nbsr_adj_pvals, gt_log_fc=log_fc$log_fc, method="NBSR")

    # use adjusted p-value to find significantly differential features.
    deseq2_ret <- evaluate(deseq2_dt$pval, gt_sig_idx, gt_not_sig_idx, thresholds)
    edgeR_ret <- evaluate(edgeR_dt$pval, gt_sig_idx, gt_not_sig_idx, thresholds)
    nbsr_ret <- evaluate(nbsr_dt$pval, gt_sig_idx, gt_not_sig_idx, thresholds)
    deseq2_ret$method <- "DESeq2"
    edgeR_ret$method <- "EdgeR"
    nbsr_ret$method <- "NBSR"
    performance_dt <- rbind(deseq2_ret, edgeR_ret, nbsr_ret)
    performance_dt$celltype <- celltype
    performance_dt$rep_no <- rep_no

    results_list <- append(results_list, list(results_dt))
    performance_list <- append(performance_list, list(performance_dt))
    
    sample_count <- dim(Y)[2]/2
    
    # Evaluate dispersion estimation.
    if (dispersion_model == "regression") {

      gt_phi <- cbind(dat_null$phi, dat_alt$phi)
      
      null_idx <- X$trt == "null"
      alt_idx <- !null_idx
      
      perturb_idxs <- log_fc$log_fc != 0
      
      avg_nbsr_phi_null <- rowMeans(nbsr_phi[,null_idx])
      sd_nbsr_phi_null <- apply(nbsr_phi[,null_idx], 1, sd)

      avg_nbsr_phi_alt <- rowMeans(nbsr_phi[,alt_idx])
      sd_nbsr_phi_alt <- apply(nbsr_phi[,alt_idx], 1, sd)
      
      nbsr_phi_null <- nbsr_phi[perturb_idxs,null_idx]
      nbsr_phi_alt <- nbsr_phi[perturb_idxs,alt_idx]

      null_phi_dt <- data.table(feature = rownames(Y)[perturb_idxs], 
                                gt = gt_phi[perturb_idxs,1],
                                edgeR_phi = edgeR_phi[perturb_idxs],
                                deseq2_phi = deseq2_phi[perturb_idxs],
                                nbsr_phi = nbsr_phi_null)
      null_phi_dt <- melt(null_phi_dt, id.vars = c("feature", "gt", "deseq2_phi", "edgeR_phi"), variable.name = "sample", value.name = "phi")
      alt_phi_dt <- data.table(feature = rownames(Y)[perturb_idxs], 
                               gt = gt_phi[perturb_idxs,2],
                               deseq2_phi = deseq2_phi[perturb_idxs],
                               edgeR_phi = edgeR_phi[perturb_idxs],
                               nbsr_phi = nbsr_phi_alt)
      alt_phi_dt <- melt(alt_phi_dt, id.vars = c("feature", "gt", "deseq2_phi", "edgeR_phi"), variable.name = "sample", value.name = "phi")
      null_phi_dt$trt <- "null"
      alt_phi_dt$trt <- "alt"
      phi_dt <- rbind(null_phi_dt, alt_phi_dt)
      feature_idx <- sample(which(perturb_idxs), 5)
      fnames <- rownames(Y)[feature_idx]
      sub_dt <- phi_dt[feature %in% fnames]
      pl <- ggplot(sub_dt, aes(trt, phi)) + theme_cowplot() +
        geom_boxplot() +
        geom_point(aes(trt, gt, shape="True phi"), size=5) +  
        geom_segment(aes(x = -Inf, xend = Inf, y = deseq2_phi, yend = deseq2_phi, color = 'DESeq2 phi'), 
                     linetype = "dashed", linewidth = 1) +  
        geom_segment(aes(x = -Inf, xend = Inf, y = edgeR_phi, yend = edgeR_phi, color = 'EdgeR phi'), 
                     linetype = "dotted", linewidth = 1) +  
        xlab("Experimental condition") + 
        ylab("Biological CV") +
        facet_grid(~ feature) +
        scale_color_manual(name = "", 
                           values = c("DESeq2 phi" = "red", "EdgeR phi" = "blue"), 
                           labels = c("DESeq2 phi" = "DESeq2 phi", "EdgeR phi" = "EdgeR phi")) +
        scale_shape_manual(name = "",
                           values = c("True phi" = 4),
                           labels = c("True phi" = "True phi"))
      pl <- pl + scale_y_log10()
      pl <- decorate_figure(pl, ylab="Dispersion (Log10)", xlab="")
      ggsave(filename = paste0(rep_path, "/biological_cv.pdf"), pl, width=12)
    }
  }
}

results <- rbindlist(results_list)
performance <- rbindlist(performance_list)

fwrite(x = results, file = paste0(data_path, "/results.csv"))
fwrite(x = performance, file = paste0(data_path, "/performance.csv"))

results <- fread(paste0(data_path, "/results.csv"))
performance <- fread(paste0(data_path, "/performance.csv"))

results$celltype2 <- gsub(pattern = "_", replacement = " ", x = stringr::str_to_title(results$celltype))
results$celltype2 <- gsub(pattern = "cd", replacement = "CD", results$celltype2)
results$celltype2 <- gsub(pattern = " CD19| CD56", replacement = "", results$celltype2)
table(results$celltype2)

performance$celltype2 <- gsub(pattern = "_", replacement = " ", x = stringr::str_to_title(performance$celltype))
performance$celltype2 <- gsub(pattern = "cd", replacement = "CD", performance$celltype2)
performance$celltype2 <- gsub(pattern = " CD19| CD56", replacement = "", performance$celltype2)
table(performance$celltype2)

# Figures will be generated in simulation_figures.R.
# figure_path <- paste0("paper/figures/swap", swap_no, "/")
# if (!dir.exists(figure_path)) {
#   dir.create(figure_path, recursive = T)
# }
# # Generate RMSE and coverage plots.
# cols <- RColorBrewer::brewer.pal(8, "Set1")
# cols2 <- RColorBrewer::brewer.pal(8, "Set2")
# pl <- ggplot(results, aes(method, rmse)) + geom_boxplot()
# pl <- pl + theme_cowplot()
# pl <- pl + ylab("RMSE") + xlab("") + theme(legend.title = element_blank())
# pl <- pl + facet_grid(~ celltype2)
# #pl <- pl + scale_fill_manual(values=cols[1:3])
# pl <- decorate_figure(pl, xlab_text = "", ylab_text = "RMSE", xtext_rotate_angle = -90)
# ggsave(filename = paste0(figure_path, "/swap", swap_no, "_rmse.pdf"), pl, width=14)
# 
# pl <- ggplot(results[method != "EdgeR"], aes(method, coverage)) + geom_boxplot()
# pl <- pl + theme_cowplot()
# pl <- pl + ylab("Coverage") + xlab("") + theme(legend.title = element_blank())
# pl <- pl + scale_fill_manual(values=cols[1:2])
# pl <- pl + geom_hline(yintercept = 0.95, linetype=2)
# pl <- pl + facet_grid(~ celltype2)
# pl <- decorate_figure(pl, xlab_text = "", ylab_text = "Coverage", xtext_rotate_angle = -90)
# ggsave(filename = paste0(figure_path, "/swap", swap_no, "_coverage.pdf"), pl, width = 14)
# 
# # Generate precision and false positive rate vs threshold.
# precision_dt <- performance[,.(avg=mean(prec), std=sd(prec)),by=.(threshold, method)]
# pl <- ggplot(precision_dt, aes(factor(threshold), avg)) + geom_point() + geom_line(aes(group=method))
# pl <- pl + theme_cowplot()
# pl <- pl + geom_errorbar(aes(ymin = avg-2*std, ymax = avg+2*std), width=0.1)
# #pl <- pl + scale_color_manual(values=cols[1:3])
# pl <- pl + ylab("Precision") + xlab("Significance threshold") + theme(legend.title = element_blank())
# pl <- pl + facet_grid(~method)
# ggsave(filename = paste0(figure_path, "/swap", swap_no, "_precision.pdf"), pl)
# 
# fpr <- performance[,.(avg=mean(fpr), std=sd(fpr)),by=.(threshold, method)]
# fnr <- performance[,.(avg=mean(fnr), std=sd(fnr)),by=.(threshold, method)]
# fpr$type <- "FPR"
# fnr$type <- "FNR"
# 
# pl <- ggplot(fpr, aes(factor(threshold), avg)) + geom_point() + geom_line(aes(group=method))
# pl <- pl + theme_cowplot()
# pl <- pl + geom_errorbar(aes(ymin = avg-2*std, ymax = avg+2*std), width=0.1)
# #pl <- pl + scale_color_manual(values=cols[1:3])
# pl <- pl + ylab("FPR") + xlab("Significance threshold") + theme(legend.title = element_blank())
# pl <- pl + facet_grid(~ method)
# ggsave(filename = paste0(figure_path, "/swap", swap_no, "_fpr.pdf"), pl)
# 
# pl <- ggplot(fnr, aes(factor(threshold), avg)) + geom_point() + geom_line(aes(group=method))
# pl <- pl + theme_cowplot(12)
# pl <- pl + geom_errorbar(aes(ymin = avg-2*std, ymax = avg+2*std), width=0.1)
# #pl <- pl + scale_color_manual(values=cols[1:3])
# pl <- pl + ylab("FNR") + xlab("Significance threshold") + theme(legend.title = element_blank())
# pl <- pl + facet_grid(~ method)
# ggsave(filename = paste0(figure_path, "/swap", swap_no, "_fnr.pdf"), pl)
# 
# # Plot FPR vs FNR as threshold is varied.
# error_rates <- rbind(fpr, fnr)
# pl <- ggplot(error_rates, aes(factor(threshold), avg, fill=type)) + geom_bar(stat = "identity", position = "stack")
# pl <- pl + facet_grid(~ method)
# pl <- pl + theme_cowplot()
# pl <- pl + scale_fill_manual(values=cols[3:4])
# pl <- pl + xlab("Significance threshold") + ylab("Error rate") + theme(legend.title = element_blank())
# ggsave(filename = paste0(figure_path, "/swap", swap_no, "_fpr_fnr.pdf"), pl)
# 
