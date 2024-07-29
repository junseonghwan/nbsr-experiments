getwd()
rm(list=ls())
library(cowplot)
library(data.table)
library(DESeq2)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(MCMCpack)
library(microRNAome)
library(patchwork)
library(seqinr)
library(stringr)
library(VennDiagram)
source("paper/functions.R")
set.seed(1)

cols <- c("NBSR+"="#00AFBB", "mRNA+"="#FC4E07") 
cols_cell_line <- c("DKO-1"="#009E73", "DKS-8"="#CC79A7", "DLD-1"="#F0E442")

UPPER_THRESHOLD <- 0.01
NS_THRESHOLD <- 0.1

FEATURES_TO_PLOT <- 8
Y_OFFSET <- -0.02
LOG2FC_THRESHOLD <- 1.5

TOP_K <- 10

M_THRESHOLD <- 4

data_path <- "data/carcinoma/"
figure_path <- paste0("paper/figures/carcinoma")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}

Y <- as.matrix(read.csv(paste0(data_path, "/Y.csv"), header = T))
X <- read.csv(paste0(data_path, "/X.csv"), header = T)
X$cell_tissue <- toupper(gsub("([a-z]+)([0-9]+)", "\\1-\\2", X$cell_tissue))
celltype1 <- "DLD-1"
celltype2 <- "DKO-1"
celltype3 <- "DKS-8"
cellline_levels <- c(celltype1, celltype2, celltype3)
X$cell_tissue <- factor(X$cell_tissue, levels = cellline_levels)

# DESeq2
se <- SummarizedExperiment(assays = list(counts = Y), colData = X)
dds2 <- DESeqDataSet(se, ~ cell_tissue)
dds2 <- DESeq(dds2, fitType = "local")
rdata <- rowData(dds2)
deseq2_mu <- assays(dds2)[["mu"]]
deseq2_phi <- rdata$dispersion
# write.table(deseq2_mu, file = paste0(data_path, "/deseq2_mu.csv"), quote=F, col.names = T, row.names=F, sep=",")
# if (WRITE_DESeq2_DISPERSION) {
#   write.table(deseq2_phi, file = paste0(data_path, "/dispersion.csv"), quote=F, col.names = F, row.names=F, sep=",")
# }
names(deseq2_phi) <- rownames(Y)

# EdgeR
d <- edgeR::DGEList(counts=Y)
d <- calcNormFactors(d)
design_mat <- model.matrix(~ cell_tissue, X)
d <- edgeR::estimateDisp(d, design = design_mat)
fit <- glmFit(d, design_mat)
edgeR_phi <- fit$dispersion
edgeR_mu <- fit$fitted.values

# NBSR
N_FEATURES <- dim(Y)[1]
R <- colSums(Y)
nbsr_pi <- read.csv(paste0(data_path, "/nbsr_pi.csv"), header = F)
nbsr_beta <- fread(paste0(data_path, "/nbsr_beta.csv"), header = F)
nbsr_beta <- matrix(nbsr_beta$V1, nrow = N_FEATURES, byrow = F)
deseq2_beta <- coef(dds2)

nbsr_phi <- as.matrix(read.csv(paste0(data_path, "/nbsr_dispersion.csv"), header = F))
eb_phi <- as.matrix(read.csv(paste0(data_path, "/eb_dispersion.csv"), header = F))

# Dispersion estimation across different methods:
plot(deseq2_phi, nbsr_phi[,1])
plot(eb_phi, nbsr_phi[,1])
plot(deseq2_phi, eb_phi)
abline(a=0, b=1)

nbsr_mu <- sweep(nbsr_pi, 2, R, "*")
rownames(nbsr_phi) <- rownames(nbsr_pi) <- rownames(nbsr_mu) <- rownames(Y)
colnames(nbsr_pi) <- colnames(nbsr_mu) <- colnames(Y)

nbsr_pi <- as.matrix(nbsr_pi)
nbsr_mu <- as.matrix(nbsr_mu)

ii <- 1
plot(nbsr_mu[,ii], Y[,ii])
points(deseq2_mu[,ii], Y[,ii], col='blue')
points(edgeR_mu[,ii], Y[,ii], col='red')
abline(a=0, b=1)

sqrt(mean((deseq2_mu - Y)^2))
sqrt(mean((edgeR_mu - Y)^2))
sqrt(mean((nbsr_mu - Y)^2)) # NBSR provides much better fit than DESeq2 and EdgeR in terms of MSE.

####### hsa-miR-10a-5p  has largest error -- how much of the MSE is explained by hsa-miR-10a-5p?

ii <- which.max(apply(abs(nbsr_mu-Y), 1, max))
print(ii)
which.max(apply(abs(deseq2_mu-Y), 1, max))
which.max(apply(abs(edgeR_mu-Y), 1, max))
plot(nbsr_mu[ii,], Y[ii,])
points(deseq2_mu[ii,], Y[ii,], col='red')
points(edgeR_mu[ii,], Y[ii,], col='blue')
abline(a = 0, b = 1)

x <- sum((deseq2_mu - Y)^2)
y <- sum((edgeR_mu - Y)^2)
z <- sum((nbsr_mu - Y)^2)
max(apply((deseq2_mu-Y)^2, 1, sum)/x)
max(apply((edgeR_mu-Y)^2, 1, sum)/y)
max(apply((nbsr_mu-Y)^2, 1, sum)/z)

generate_results_table <- function(celltype1, celltype2)
{
  celltype1_short <- celltype1 <- as.character(celltype1)
  celltype2_short <- celltype2 <- as.character(celltype2)
  #celltypes <- c(celltype1, celltype2)
  contrast_path <- paste0(data_path, "/", celltype1, "_", celltype2)
  if (!dir.exists(contrast_path)) {
    dir.create(contrast_path, recursive = T)
  }
  
  celltype1_idx <- which(X$cell_tissue == celltype1)
  celltype2_idx <- which(X$cell_tissue == celltype2)
  
  res2 <- results(dds2, contrast = c("cell_tissue", celltype1, celltype2))
  
  vec <- rep(0, length(colnames(design_mat)))
  names(vec) <- colnames(design_mat)
  colname1 <- paste0("cell_tissue", celltype1)
  colname2 <- paste0("cell_tissue", celltype2)
  if (colname1 %in% names(vec)) {
    vec[colname1] <- 1
  }
  if (colname2 %in% names(vec)) {
    vec[colname2] <- -1
  }
  lrt <- glmLRT(fit, contrast=vec)
  edgeR_results <- topTags(lrt, n=Inf, sort.by = "none")
  
  nbsr_log2_fc <- log2(nbsr_pi[,celltype1_idx[1]]) - log2(nbsr_pi[,celltype2_idx[1]])
  nbsr_logRR <- as.matrix(read.csv(paste0(contrast_path, "/nbsr_logRR.csv"), header = F))
  nbsr_logRR_sd <- as.matrix(read.csv(paste0(contrast_path, "/nbsr_logRR_sd.csv"), header = F))
  nbsr_stats <- nbsr_logRR[1,] / nbsr_logRR_sd[1,]
  nbsr_pvals <- 2*pnorm(abs(nbsr_stats), lower.tail = FALSE)
  nbsr_adj_pvals <- p.adjust(nbsr_pvals, method = "BH")
  nbsr_res <- data.table(features=rownames(Y), 
                         log2FoldChange = nbsr_log2_fc, 
                         padj=nbsr_adj_pvals,
                         pi1 = nbsr_pi[,which(X$cell_tissue == celltype1)[1]],
                         pi2 = nbsr_pi[,which(X$cell_tissue == celltype2)[1]])
  
  write.csv(res2, paste0(contrast_path, "/DESeq2_res.csv"), quote = F, row.names = T)
  write.csv(edgeR_results, paste0(contrast_path, "/edgeR_res.csv"), quote = F, row.names = T)
  write.csv(nbsr_res, paste0(contrast_path, "/NBSR_res.csv"), quote = F, row.names = T)
  
  res2_dt <- data.table(features = rownames(res2), 
                        log2FoldChange=res2$log2FoldChange,
                        padj=res2$padj)
  edgeR_dt <- data.table(features = rownames(edgeR_results$table), 
                         log2FoldChange=edgeR_results$table$logFC,
                         padj=edgeR_results$table$FDR)
  return(list("DESeq2"=res2_dt, "EdgeR"=edgeR_dt, "NBSR"=nbsr_res))
}
annotate_conditions <- function(res, cond_A, cond_B, factor_levels)
{
  res[,condition_A := factor(cond_A, levels = factor_levels)]
  res[,condition_B := factor(cond_B, levels = factor_levels)]
}

results_list_12 <- generate_results_table(celltype1 = celltype1, celltype2 = celltype2)
results_list_13 <- generate_results_table(celltype1 = celltype1, celltype2 = celltype3)
results_list_23 <- generate_results_table(celltype1 = celltype2, celltype2 = celltype3)
results_list_12 <- lapply(results_list_12, annotate_conditions, cond_A=celltype1, cond_B=celltype2, factor_levels=cellline_levels)
results_list_23 <- lapply(results_list_23, annotate_conditions, cond_A=celltype2, cond_B=celltype3, factor_levels=cellline_levels)
results_list_13 <- lapply(results_list_13, annotate_conditions, cond_A=celltype1, cond_B=celltype3, factor_levels=cellline_levels)

nbsr_res_12 <- results_list_12[["NBSR"]]
nbsr_res_13 <- results_list_13[["NBSR"]]
nbsr_res_23 <- results_list_23[["NBSR"]]

# Composition plot for each cell line.
celltype1_idx <- which(X$cell_tissue == celltype1)
celltype2_idx <- which(X$cell_tissue == celltype2)
celltype3_idx <- which(X$cell_tissue == celltype3)
celltype1_short <- celltype1
celltype2_short <- celltype2
celltype3_short <- celltype3

pi_A <- nbsr_pi[,celltype1_idx[1]]
pi_B <- nbsr_pi[,celltype2_idx[1]]
pi_C <- nbsr_pi[,celltype3_idx[1]]

head(pi_A)
head(pi_B)
head(pi_C)
cumsum(sort(pi_A, decreasing = T))[1:TOP_K]
cumsum(sort(pi_B, decreasing = T))[1:TOP_K]
cumsum(sort(pi_C, decreasing = T))[1:TOP_K]

# 1. Composition bar plot to compare samples.
# Convert to data frames and add a condition label
dt_A <- data.table(features = names(pi_A), probability = pi_A, condition = celltype1_short)
dt_B <- data.table(features = names(pi_B), probability = pi_B, condition = celltype2_short)
dt_C <- data.table(features = names(pi_C), probability = pi_C, condition = celltype3_short)

dt_top_K_A <- dt_A %>% top_n(TOP_K, probability) %>% arrange(desc(probability))
dt_top_K_B <- dt_B %>% top_n(TOP_K, probability) %>% arrange(desc(probability))
dt_top_K_C <- dt_C %>% top_n(TOP_K, probability) %>% arrange(desc(probability))
unique_features <- unique(c(dt_top_K_A$feature, dt_top_K_B$feature, dt_top_K_C$feature))

dt_combined <- bind_rows(dt_A[features %in% unique_features],
                         dt_B[features %in% unique_features],
                         dt_C[features %in% unique_features])
dt_combined$miRNA_short <- gsub("hsa-", "", dt_combined$features)
# Order the miRNAs by probability in A.
temp_A <- dt_combined[condition == celltype1_short]
temp_A <- temp_A[order(probability, decreasing = T)]
dt_combined$miRNA_short <- factor(dt_combined$miRNA_short, levels = temp_A$miRNA_short)
dt_combined[,condition := factor(condition, levels = cellline_levels)]

unique(dt_combined$features)
# Create a grouped bar chart
pl <- ggplot(dt_combined, aes(x = miRNA_short, y = probability)) +
  geom_bar(aes(fill = condition), stat = "identity", position = position_dodge(), width = 0.9) +
  labs(x = "", y = "Estimated proportion") +
  theme_minimal() +
  guides(fill = guide_legend(title=NULL))
pl <- pl + scale_fill_manual(values = cols_cell_line)
#filename <- paste0(figure_path, "/", TOP_K, "_miRNAs_log2fc.pdf")
#ggsave(filename = filename, pl, width=18, height = 6)

# Now, annotate p-values.
dt_AB <- left_join(dt_A, nbsr_res_12)
dt_BC <- left_join(dt_B, nbsr_res_23)
dt_AC <- left_join(dt_C, nbsr_res_13)

joined_dt <- bind_rows(dt_AB, dt_BC, dt_AC)
joined_dt <- subset(joined_dt, joined_dt$features %in% unique_features)
joined_dt$miRNA_short <- gsub("hsa-", "", joined_dt$features)
joined_dt$miRNA_short <- factor(joined_dt$miRNA_short, levels = temp_A$miRNA_short)
joined_dt[,condition := factor(condition, levels = cellline_levels)]
joined_dt[,condition_A := factor(condition_A, levels = cellline_levels)]
joined_dt[,condition_B := factor(condition_B, levels = cellline_levels)]
joined_dt$x <- as.numeric(joined_dt$miRNA_short)
joined_dt$y <- as.numeric(joined_dt$probability)

# Get max y for each feature
max_y_dt <- joined_dt[,.(max_y=max(y)),by=.(features)]

# Generate stats to plot for each pair of condition.
get_pval_category <- function(pval) {
  ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01, "**",
                ifelse(pval < 0.05, "*", "NS")))
}
annotate_plot <- function(pl_to_decorate, celltype, x_offset1, x_offset2, y_offset)
{
  top_K_estimated_prop <- joined_dt[condition == celltype,.(x, y=max(y), log2FoldChange, padj, condition, condition_A, condition_B),by=.(features)]
  top_K_estimated_prop$xbegin <- top_K_estimated_prop$x + x_offset1
  top_K_estimated_prop$xend <- top_K_estimated_prop$x + x_offset2
  #top_K_estimated_prop$ybegin <- top_K_estimated_prop$yend <- top_K_estimated_prop$y
  top_K_estimated_prop$ybegin <- 0
  top_K_estimated_prop <- left_join(top_K_estimated_prop, max_y_dt)
  top_K_estimated_prop$yend <- top_K_estimated_prop$max_y + y_offset
  
  top_K_estimated_prop$pval_cat <- get_pval_category(top_K_estimated_prop$padj)
  top_K_estimated_prop$pval_cat <- factor(top_K_estimated_prop$pval_cat, levels = c("***", "**", "*", "NS"))
  top_K_estimated_prop$significance <- abs(top_K_estimated_prop$log2FoldChange) > LOG2FC_THRESHOLD
  top_K_estimated_prop <- top_K_estimated_prop[pval_cat != "NS"]
  
  pl_to_decorate <- pl_to_decorate + geom_segment(data=top_K_estimated_prop, aes(x=xbegin, xend=xend, y=yend + 0.01, yend=yend + 0.01), inherit.aes = FALSE)  
  pl_to_decorate <- pl_to_decorate + geom_text(data=top_K_estimated_prop, aes(x=(xbegin+xend)/2, y=yend + 0.015, label=pval_cat), 
                                               col = "black",
                                               inherit.aes = FALSE)
  pl_to_decorate <- pl_to_decorate + theme(plot.title = element_text(size = 20))  # Enlarge title size
  return(pl_to_decorate)
}

pl2 <- annotate_plot(pl, celltype1, -0.45, 0.15, 0)
pl3 <- annotate_plot(pl2, celltype2, -0.15, 0.5, 0.025)
pl4 <- annotate_plot(pl3, celltype3, -0.45, 0.5, 0.05)
pl4 <- pl4 + ggtitle(paste0(TOP_K, " miRNAs with highest estimated proportion in either cell type"))
pl4 <- pl4 + theme(plot.title = element_text(size = 20))  # Enlarge title size
pl4 <- pl4 + coord_cartesian(clip="off")
pl4 <- decorate_figure(pl4, xlab_text = "", ylab_text = "Estimated Proportion", xtext_rotate_angle = -45)
filename <- paste0(figure_path, "/", TOP_K, "_miRNAs_pvals.pdf")
ggsave(filename = filename, pl4, width=9, height = 6)

# DLD-1 : DKO-1
fname <- "hsa-miR-10a-5p"
#fname <- "hsa-miR-200b-3p"
#fname <- "hsa-let-7i-5p"

nbsr_pi[fname, celltype1_idx[1]]
nbsr_pi[fname, celltype2_idx[1]]
nbsr_pi[fname, celltype3_idx[1]]

results_list_12$DESeq2[features == fname,]
results_list_12$EdgeR[features == fname,]
results_list_12$NBSR[features == fname,]

# DKO-1 : DKS-8
results_list_23$DESeq2[features == fname,]
results_list_23$EdgeR[features == fname,]
results_list_23$NBSR[features == fname,]

# DLD-1 : DKS-8
results_list_13$DESeq2[features == fname,]
results_list_13$EdgeR[features == fname,]
results_list_13$NBSR[features == fname,]

dt_combined[,sum(probability),by=.(condition)]

#######

# 2. Examine miRNAs proportions and normalized counts across three groups.

#######

table(X$organ, X$cell_tissue)

# Normalized counts showing the differing spread across the conditions.
deseq2_normalized_cts <- counts(dds2, normalized=TRUE)
X_dt <- as.data.table(X[,c("sample_id", "cell_tissue")])
X_dt$short_name <- X_dt$cell_tissue

# Get normalized counts.
cts_dt <- data.table(deseq2_normalized_cts,
                     features=rownames(Y))
cts_dt_melt <- melt(cts_dt, id.vars = "features", variable.name = "sample_id", value.name = "value")
cts_dt_melt_join <- dplyr::left_join(cts_dt_melt, X_dt, by = "sample_id")

# Props
props <- apply(Y, 2, normalize)
props_dt <- data.table(t(props), sample_id=X$sample_id)
props_dt_melt <- melt(props_dt, id.vars = c("sample_id"), variable.name = "features")
props_dt_melt_join <- dplyr::left_join(props_dt_melt, X_dt, by = "sample_id")

# Uses objects in the global environment.
# Plot the proportions and NBSR adjusted p-values along with
# normalized counts and DESeq2 adjusted p-values to compare all three groups.
generate_prop_norm_cts_plots <- function(f_idxs, fig_output_path, file_prefix="", file_suffix="")
{
  fnames <- rownames(Y)[f_idxs]
  
  for (fname in fnames)
  {
    temp_cts_dt <- cts_dt_melt_join[features == fname]
    temp_props_dt <- props_dt_melt_join[features == fname]
    
    estimated_prop <- nbsr_pi[rownames(nbsr_pi) == fname,
                              c(celltype1_idx[1], 
                                celltype2_idx[1],
                                celltype3_idx[1])]
    temp_pi_dt <- data.table(condition=cellline_levels, 
                             pi=estimated_prop)
    temp_pi_avg_dt <- left_join(temp_pi_dt, temp_cts_dt[,.(avg_norm_cts=mean(value),condition=cell_tissue),by=.(cell_tissue)])
    
    # p-values for DESeq2
    deseq2_pvals_dt <- bind_rows(results_list_12$DESeq2[features == fname],
                                 results_list_23$DESeq2[features == fname],
                                 results_list_13$DESeq2[features == fname])
    deseq2_pvals_dt$xbegin <- as.numeric(deseq2_pvals_dt$condition_A)
    deseq2_pvals_dt$xend <- as.numeric(deseq2_pvals_dt$condition_B)
    deseq2_pvals_dt$yend <- deseq2_pvals_dt$ybegin <- max(temp_cts_dt$value) * (1 + 1:length(cellline_levels)/8)
    
    # p-values for NBSR
    nbsr_pvals_dt <- bind_rows(results_list_12$NBSR[features == fname],
                               results_list_23$NBSR[features == fname],
                               results_list_13$NBSR[features == fname])
    nbsr_pvals_dt$xbegin <- as.numeric(nbsr_pvals_dt$condition_A)
    nbsr_pvals_dt$xend <- as.numeric(nbsr_pvals_dt$condition_B)
    nbsr_pvals_dt$yend <- nbsr_pvals_dt$ybegin <- max(temp_props_dt$value) * (1 + 1:length(cellline_levels)/8)
    
    # pl_norm_cts <- ggboxplot(data = temp_cts_dt,
    #                          x = "cell_tissue", y = "value",
    #                          add="jitter")
    pl_norm_cts <- ggplot(data = temp_cts_dt, aes(cell_tissue, y = value)) + 
      geom_boxplot() + geom_jitter() + theme_minimal()
    pl_norm_cts <- pl_norm_cts + geom_segment(data = deseq2_pvals_dt,
                                              aes(x=xbegin, xend=xend, 
                                                  y=ybegin, yend=yend))
    pl_norm_cts <- pl_norm_cts + geom_text(data = deseq2_pvals_dt,
                                           aes(x=(xbegin+xend)/2, 
                                               y=ybegin*1.05,
                                               label=signif(padj, 3)))
    pl_norm_cts <- pl_norm_cts + ggtitle("DESeq2 adj. p-values")
    pl_norm_cts <- pl_norm_cts + geom_point(data = temp_pi_avg_dt,
                                            aes(condition, avg_norm_cts, col="Avg of normalized cts")) +
      scale_color_manual(values=c("Avg of normalized cts"="red")) +
      theme(legend.title = element_blank(), legend.position = "bottom",
            legend.text = element_text(size = 16),
            legend.box.margin = margin(-20, 0, 0, 0))
    pl_norm_cts <- decorate_figure(pl_norm_cts, xlab="", ylab_text = "Normalized counts")
    
    # pl_props <- ggboxplot(data = temp_props_dt,
    #                       x = "cell_tissue", y = "value",
    #                       add="jitter")
    pl_props <- ggplot(data = temp_props_dt, aes(cell_tissue, y = value)) + 
      geom_boxplot() + geom_jitter() + theme_minimal()
    pl_props <- pl_props + geom_segment(data = nbsr_pvals_dt,
                                        aes(x=xbegin, xend=xend, 
                                            y=ybegin, yend=yend))
    pl_props <- pl_props + geom_text(data = nbsr_pvals_dt,
                                     aes(x=(xbegin+xend)/2, 
                                         y=ybegin*1.05,
                                         label=signif(padj, 3)))
    pl_props <- pl_props + geom_point(data = temp_pi_avg_dt,
                                      aes(condition, pi, col="NBSR estimated proportion")) +
      scale_color_manual(values=c("NBSR estimated proportion"="red")) +
      theme(legend.title = element_blank(), legend.position = "bottom",
            legend.text = element_text(size = 16),
            legend.box.margin = margin(-20, 0, 0, 0))
    pl_props <- pl_props + ggtitle("NBSR adj. p-values")
    pl_props <- decorate_figure(pl_props, xlab="", ylab_text = "Sample proportion")
    
    fname_short <- gsub("hsa-", "", fname)
    pl <- pl_props + pl_norm_cts
    pl <- pl + plot_annotation(title = fname_short, theme = theme(plot.title=element_text(size=24, hjust=0.5)))
    filename <- paste0(file_prefix, "prop_norm_cts_", gsub("/", "|", fname_short), file_suffix, ".pdf")
    print(filename)
    ggsave(filename = paste0(fig_output_path, "/", filename), pl, width = 8, height = 5)
  }
}

fname <- "hsa-miR-10a-5p"
f_idxs <- which(grepl(fname, rownames(Y)))
generate_prop_norm_cts_plots(f_idxs, figure_path)

results_list_12$DESeq2[features == fname]
results_list_12$EdgeR[features == fname]
results_list_12$NBSR[features == fname]

f_idxs <- which(rownames(Y) %in% unique_features)
generate_prop_norm_cts_plots(f_idxs, figure_path, file_prefix = paste0("Top", TOP_K, "_"))

# 3. MA plots: pairwise comparisons.

avg_normalized_cts <- cts_dt_melt_join[,.(M_norm_cts=mean(value)),by=.(features)]
avg_normalized_cts[,short_name := gsub("hsa-", "", features)]

avg_props <- props_dt_melt_join[,.(M_props=mean(value)),by=.(features)]
avg_props[,short_name := gsub("hsa-", "", features)]

avg_norm_cts_props <- left_join(avg_normalized_cts, avg_props)

# Retrieve log2fc and FDRs for all three methods.
generate_figures <- function(res, file_suffix="", nudge_x=0, nudge_y=0)
{
  contrast_path <- paste0(figure_path, "/", res$NBSR$condition_A[1], "_", res$NBSR$condition_B[1])
  if (!dir.exists(contrast_path)) {
    dir.create(contrast_path, recursive = T)
  }
  subfigure_path <- paste0(contrast_path, "/subfigures")
  if (!dir.exists(subfigure_path)) {
    dir.create(subfigure_path, recursive = T)
  }
  
  log2fc_dt <- data.table(features=rownames(Y), 
                          NBSR = res$NBSR$log2FoldChange, 
                          DESeq2 = res$DESeq2$log2FoldChange, 
                          EdgeR = res$EdgeR$log2FoldChange)
  adj_pval_dt <- data.table(features=rownames(Y), 
                            NBSR_pval = res$NBSR$padj, 
                            DESeq2_pval = res$DESeq2$padj, 
                            EdgeR_pval = res$EdgeR$padj)
  
  avg_norm_cts_props_log2fc <- left_join(avg_norm_cts_props, log2fc_dt)
  
  # dim(adj_pval_dt[NBSR_pval <= UPPER_THRESHOLD & DESeq2_pval <= UPPER_THRESHOLD])
  # dim(adj_pval_dt[NBSR_pval <= UPPER_THRESHOLD & DESeq2_pval > NS_THRESHOLD])
  # 
  # dim(adj_pval_dt[NBSR_pval <= UPPER_THRESHOLD & EdgeR_pval <= UPPER_THRESHOLD])
  # dim(adj_pval_dt[NBSR_pval <= UPPER_THRESHOLD & EdgeR_pval > NS_THRESHOLD])
  
  # Find all miRNAs where the methods differ significantly.
  nbsr_adj_pval_dt <- adj_pval_dt[NBSR_pval <= UPPER_THRESHOLD & DESeq2_pval > NS_THRESHOLD & EdgeR_pval > NS_THRESHOLD,]
  mRNA_adj_pval_dt <- adj_pval_dt[NBSR_pval > NS_THRESHOLD & (DESeq2_pval <= UPPER_THRESHOLD | EdgeR_pval <= UPPER_THRESHOLD),]
  nbsr_adj_pval_dt[,method := "NBSR+"]
  mRNA_adj_pval_dt[,method := "mRNA+"]
  
  sub_adj_pval_dt <- bind_rows(nbsr_adj_pval_dt, mRNA_adj_pval_dt)
  sub_log2fc_dt <- log2fc_dt[features %in% sub_adj_pval_dt$features]
  sub_log2fc_dt <- left_join(sub_log2fc_dt, sub_adj_pval_dt)
  nbsr_mRNA_diff_plot_dat <- left_join(sub_log2fc_dt, avg_norm_cts_props)
  
  # NBSR vs DESeq2.
  sub_pl_dat_DESeq2 <- avg_norm_cts_props_log2fc[abs(NBSR - DESeq2) > 3]
  # Log2FC plot.
  pl <- ggplot(avg_norm_cts_props_log2fc, aes(NBSR, DESeq2)) + 
    geom_point(alpha=0.1) + theme_minimal()
  pl <- pl + ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B))
  #   labs(subtitle = "Comparison of Log2 FC")
  pl <- pl + geom_abline(slope = 1, intercept = 0)
  pl <- pl + geom_point(data = sub_pl_dat_DESeq2, aes(NBSR, DESeq2), col="red", alpha=0.5)
  pl <- pl + geom_text_repel(data = sub_pl_dat_DESeq2[log1p(M_norm_cts) > M_THRESHOLD],
                             aes(NBSR, DESeq2, label = short_name), nudge_y = 1)
  pl <- decorate_figure(pl, xlab_text = "NBSR Log2FC", ylab_text = "DESeq2 Log2FC")
  ggsave(filename = paste0(contrast_path, "/log2FC_NBSR_DESeq2", file_suffix, ".pdf"), pl)
  
  # MA-plots for NBSR.
  ma_pl <- ggplot(avg_norm_cts_props_log2fc, aes(log1p(M_norm_cts), NBSR)) + 
    geom_point(alpha=0.1) +
    theme_minimal()
  ma_pl <- ma_pl + geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", linewidth = 0.5, col ='red')
  ma_pl <- ma_pl + ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B))
  ma_pl1 <- ma_pl + geom_point(data=sub_pl_dat_DESeq2, aes(log1p(M_norm_cts), NBSR), col="red", alpha=0.5)
  ma_pl1 <- ma_pl1 + geom_text_repel(data=sub_pl_dat_DESeq2[log1p(M_norm_cts) > 4],
                                     aes(log1p(M_norm_cts), NBSR, label=short_name), 
                                     size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl1 <- decorate_figure(ma_pl1, xlab_text = "Log1p(M)", 
                            ylab_text = "NBSR Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_NBSR_DESeq2", file_suffix, ".pdf"), ma_pl1)
  
  ma_pl2 <- ma_pl + geom_point(data = nbsr_mRNA_diff_plot_dat, 
                               aes(log1p(M_norm_cts), NBSR, col=method)) +
    scale_color_manual(values = cols) +
    theme(legend.title = element_blank(), legend.position = c(0.85, 0.1)) +
    ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B)) +
    geom_text_repel(data=nbsr_mRNA_diff_plot_dat[log1p(M_norm_cts) > 4 & method == "mRNA+"],
                    aes(log1p(M_norm_cts), NBSR, label=short_name),
                    size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl2 <- decorate_figure(ma_pl2, xlab_text = "Log1p(M)", 
                            ylab_text = "NBSR Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_NBSR_diff", file_suffix, ".pdf"), ma_pl2)
  
  # MA-plots for DESeq2.
  ma_pl <- ggplot(avg_norm_cts_props_log2fc, aes(log1p(M_norm_cts), DESeq2)) + 
    geom_point(alpha=0.1) + 
    theme_minimal() 
  ma_pl <- ma_pl + geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", linewidth = 0.5, col ='red')
  ma_pl <- ma_pl + ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B))
  
  ma_pl1 <- ma_pl + geom_point(data=sub_pl_dat_DESeq2, aes(log1p(M_norm_cts), DESeq2), col="red", alpha=0.5)
  ma_pl1 <- ma_pl1 + geom_text_repel(data=sub_pl_dat_DESeq2[log1p(M_norm_cts) > 4],
                                     aes(log1p(M_norm_cts), DESeq2, label=short_name),
                                     size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl1 <- decorate_figure(ma_pl1, xlab_text = "Log1p(M)", 
                            ylab_text = "DESeq2 Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_DESeq2", file_suffix, ".pdf"), ma_pl1)
  
  ma_pl2 <- ma_pl + geom_point(data = nbsr_mRNA_diff_plot_dat, 
                               aes(log1p(M_norm_cts), DESeq2, col=method)) +
    scale_color_manual(values = cols) +
    theme(legend.title = element_blank(), legend.position = c(0.85, 0.1)) +
    ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B)) + 
    geom_text_repel(data=nbsr_mRNA_diff_plot_dat[log1p(M_norm_cts) > 4 & method == "mRNA+"],
                    aes(log1p(M_norm_cts), DESeq2, label=short_name),
                    size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl2 <- decorate_figure(ma_pl2, xlab_text = "Log1p(M)", 
                            ylab_text = "DESeq2 Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_DESeq2_diff", file_suffix, ".pdf"), ma_pl2)
  
  # NBSR vs EdgeR
  sub_pl_dat_edgeR <- avg_norm_cts_props_log2fc[abs(NBSR - EdgeR) > 3]
  # Log2FC plot.
  pl <- ggplot(avg_norm_cts_props_log2fc, aes(NBSR, EdgeR)) + 
    geom_point(alpha=0.1) + theme_minimal()
  pl <- pl + geom_abline(slope = 1, intercept = 0)
  pl <- pl + ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B))
  #  labs(subtitle = "Comparison of Log2 FC")
  pl <- pl + geom_point(data = sub_pl_dat_edgeR, aes(NBSR, EdgeR), col="red", alpha=0.5)
  pl <- pl + geom_text_repel(data = sub_pl_dat_edgeR[log1p(M_norm_cts) > M_THRESHOLD],
                             aes(NBSR, EdgeR, label = short_name), nudge_y = 1)
  pl <- decorate_figure(pl, xlab_text = "NBSR Log2FC", ylab_text = "EdgeR Log2FC")
  ggsave(filename = paste0(contrast_path, "/log2FC_NBSR_EdgeR", file_suffix, ".pdf"), pl)
  
  # MA plot for NBSR.
  ma_pl <- ggplot(avg_norm_cts_props_log2fc, aes(log1p(M_norm_cts), NBSR)) + 
    geom_point(alpha=0.1) +
    theme_minimal()
  ma_pl <- ma_pl + geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", linewidth = 0.5, col ='red')
  ma_pl <- ma_pl + ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B))
  
  ma_pl1 <- ma_pl + geom_point(data=sub_pl_dat_edgeR, aes(log1p(M_norm_cts), NBSR), col="red", alpha=0.5)
  ma_pl1 <- ma_pl1 + geom_text_repel(data=sub_pl_dat_edgeR[log1p(M_norm_cts) > 4],
                                     aes(log1p(M_norm_cts), NBSR, label=short_name),
                                     size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl1 <- decorate_figure(ma_pl1, xlab_text = "Log1p(M)", 
                            ylab_text = "EdgeR Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_NBSR_EdgeR", file_suffix, ".pdf"), ma_pl1)
  
  # MA-plots for EdgeR
  ma_pl <- ggplot(avg_norm_cts_props_log2fc, aes(log1p(M_norm_cts), EdgeR)) + 
    geom_point(alpha=0.1) + 
    theme_minimal() 
  ma_pl <- ma_pl + geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", linewidth = 0.5, col ='red')
  
  ma_pl1 <- ma_pl + geom_point(data=sub_pl_dat_edgeR, aes(log1p(M_norm_cts), EdgeR), col="red", alpha=0.5)
  ma_pl1 <- ma_pl1 + geom_text_repel(data=sub_pl_dat_edgeR[log1p(M_norm_cts) > 4],
                                     aes(log1p(M_norm_cts), EdgeR, label=short_name),
                                     size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl1 <- decorate_figure(ma_pl1, xlab_text = "Log1p(M)", 
                            ylab_text = "NBSR Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_EdgeR", file_suffix, ".pdf"), ma_pl1)
  
  ma_pl2 <- ma_pl + geom_point(data = nbsr_mRNA_diff_plot_dat, 
                               aes(log1p(M_norm_cts), EdgeR, col=method)) +
    scale_color_manual(values = cols) +
    theme(legend.title = element_blank(), legend.position = c(0.85, 0.1)) +
    ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B)) + 
    geom_text_repel(data=nbsr_mRNA_diff_plot_dat[log1p(M_norm_cts) > 4 & method == "mRNA+"],
                    aes(log1p(M_norm_cts), EdgeR, label=short_name),
                    size=5, nudge_x = nudge_x, nudge_y = nudge_y)
  ma_pl2 <- decorate_figure(ma_pl2, xlab_text = "Log1p(M)", 
                            ylab_text = "EdgeR Log2FC")
  ggsave(filename = paste0(contrast_path, "/MA_plot_EdgeR_diff", file_suffix, ".pdf"), ma_pl2)
  
  # Log2FC plot for DESeq2 vs EdgeR.
  pl <- ggplot(avg_norm_cts_props_log2fc, aes(EdgeR, DESeq2)) + 
    geom_point(alpha=0.1) + theme_minimal()
  pl <- pl + geom_abline(slope = 1, intercept = 0)
  pl <- pl + ggtitle(paste0(res$NBSR$condition_A, " vs ", res$NBSR$condition_B))
  #labs(subtitle = "Comparison of Log2 FC")
  pl <- decorate_figure(pl, xlab_text = "EdgeR", ylab_text = "DESeq2")
  ggsave(filename = paste0(contrast_path, "/log2FC_DESeq2_EdgeR", file_suffix, ".pdf"), pl)
}
generate_figures(results_list_12)
generate_figures(results_list_23)
generate_figures(results_list_13)


# Generate normalized counts and sample proportion plot for points of interest.
get_contrast_path <- function(res) 
{
  contrast_path <- paste0(figure_path, "/", res$NBSR$condition_A[1], "_", res$NBSR$condition_B[1])
  return(contrast_path)
}
# 1 vs 2
fnames <- c("hsa-miR-708-3p", "hsa-miR-4723-5p")
f_idxs <- which(rownames(Y) %in% fnames)
contrast_path <- get_contrast_path(results_list_12)
contrast_subfig_path <- paste0(contrast_path, "/subfigures")
generate_prop_norm_cts_plots(f_idxs, contrast_subfig_path)

deseq2_normalized_cts[fnames[1],celltype1_idx]
deseq2_normalized_cts[fnames[1],celltype2_idx]

# 1 vs 3
fnames <- c("hsa-miR-520d-3p")
f_idxs <- which(rownames(Y) %in% fnames)
contrast_path <- get_contrast_path(results_list_13)
contrast_subfig_path <- paste0(contrast_path, "/subfigures")
generate_prop_norm_cts_plots(f_idxs, contrast_subfig_path)

# 2 vs 3
fnames <- c("hsa-miR-708-3p", "hsa-miR-520d-3p")
f_idxs <- which(rownames(Y) %in% fnames)
contrast_path <- get_contrast_path(results_list_23)
contrast_subfig_path <- paste0(contrast_path, "/subfigures")
generate_prop_norm_cts_plots(f_idxs, contrast_subfig_path)
