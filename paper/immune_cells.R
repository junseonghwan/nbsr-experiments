rm(list=ls())
library(cowplot)
library(data.table)
library(DESeq2)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(MCMCpack)
library(microRNAome)
library(patchwork)
library(seqinr)
library(stringr)
library(VennDiagram)
source("paper/functions.R")
set.seed(1)

datapath <- "data/immune/"
if (!dir.exists(datapath)) {
  dir.create(datapath, recursive = T)
}
figure_path <- paste0("paper/figures/immune/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}
subfigure_path <- paste0(figure_path, "/subfigures")
if (!dir.exists(subfigure_path)) {
  dir.create(subfigure_path, recursive = T)
}

# Load data from scratch and apply filtering.
data("microRNAome")

### microRNAome
table(microRNAome$cell_tissue)
table(microRNAome$study_id)
temp <- microRNAome[,microRNAome$study_id == 89]
sort(table(temp$cell_tissue))
table(temp$sequencer)
table(temp$cell_tissue, temp$organ)
celltype1 <- "T_lymphocyte_CD8"
celltype2 <- "B_lymphocyte_CD19"
sample_i <- grepl(pattern = paste0("(", celltype1, "|", celltype2, ")"), temp$cell_tissue)
sum(sample_i)
table(temp$cell_tissue[sample_i])
table(temp$organ[sample_i])
dat <- temp[,sample_i]
i_rm <- (colSums(assay(dat)) < 1e6)
dat <- dat[,!i_rm]
dim(dat)
cts_pre_filtering <- colSums(assay(dat))

celltype1_short <- gsub("_lymphocyte_", " ", celltype1)
celltype2_short <- gsub("_lymphocyte_", " ", celltype2)

cts <- assay(dat)
#Y <- cts[!(rowMeans(cts == 0) == 1),]
FEATURE_EXPRESSION_MIN_SAMPLES <- 3
Y <- cts[(rowSums(cts > 0) >= FEATURE_EXPRESSION_MIN_SAMPLES),]
X <- colData(dat)
dim(Y)
dim(X)

table(X$cell_tissue)

cts_pre_filtering - colSums(Y)

# Uncomment below to overwrite the data.
#write.table(X, file = paste0(out_path, "/X.csv"), quote = F, col.names = T, row.names = F, sep=",")
#write.table(Y, file = paste0(out_path, "/Y.csv"), quote = F, col.names = T, row.names = T, sep=",")

celltype1_idx <- which(X$cell_tissue == celltype1)
celltype2_idx <- which(X$cell_tissue == celltype2)

table(X$cell_tissue)
se <- SummarizedExperiment(assays = list(counts = Y), colData = X)
dds2 <- DESeqDataSet(se, ~ cell_tissue)
dds2 <- DESeq(dds2, fitType = "local")
res2 <- results(dds2, contrast = c("cell_tissue", celltype1, celltype2))
rdata <- rowData(dds2)
deseq2_mu <- assays(dds2)[["mu"]]
deseq2_phi <- rdata$dispersion
write.table(deseq2_mu, file = paste0(datapath, "/deseq2_mu.csv"), quote=F, col.names = T, row.names=F, sep=",")
names(deseq2_phi) <- rownames(Y)

d <- edgeR::DGEList(counts=Y)
d <- calcNormFactors(d)
design_mat <- model.matrix(~ cell_tissue, X)
d <- edgeR::estimateDisp(d, design = design_mat)
fit <- glmFit(d, design_mat)
lrt <- glmLRT(fit, contrast=c(0, 1))
edgeR_results <- topTags(lrt, n=Inf, sort.by = "none")
edgeR_phi <- fit$dispersion
edgeR_mu <- fit$fitted.values

R <- colSums(Y)
nbsr_pi <- read.csv(paste0(out_path, "/nbsr_pi.csv"), header = F)
nbsr_phi <- as.matrix(read.csv(paste0(out_path, "/nbsr_dispersion.csv"), header = F))
dim(nbsr_pi)
nbsr_mu <- sweep(nbsr_pi, 2, R, "*")
plot(nbsr_mu[,1], Y[,1])
points(deseq2_mu[,1], Y[,1], col='red')
points(edgeR_mu[,1], Y[,1], col='blue')

rownames(nbsr_phi) <- rownames(nbsr_pi) <- rownames(nbsr_mu) <- rownames(Y)
colnames(nbsr_phi) <- colnames(nbsr_pi) <- colnames(nbsr_mu) <- colnames(Y)

nbsr_pi <- as.matrix(nbsr_pi)
nbsr_mu <- as.matrix(nbsr_mu)
nbsr_log2_fc <- log2(nbsr_pi[,which(X$cell_tissue == celltype1)[1]]) - log2(nbsr_pi[,which(X$cell_tissue == celltype2)[1]])
results_path <- paste0(out_path, paste0(celltype1, "_", celltype2))
nbsr_logRR <- as.matrix(read.csv(paste0(results_path, "/nbsr_logRR.csv"), header = F))
nbsr_logRR_sd <- as.matrix(read.csv(paste0(results_path, "/nbsr_logRR_sd.csv"), header = F))
nbsr_stats <- nbsr_logRR[1,] / nbsr_logRR_sd[1,]
nbsr_pvals <- 2*pnorm(abs(nbsr_stats), lower.tail = FALSE)
nbsr_adj_pvals <- p.adjust(nbsr_pvals, method = "BH")
nbsr_res <- data.table(features=rownames(Y), log2FoldChange = nbsr_log2_fc, padj=nbsr_adj_pvals)

sqrt(mean((deseq2_mu - Y)^2))
sqrt(mean((edgeR_mu - Y)^2))
sqrt(mean((nbsr_mu - Y)^2))

sum(res2$padj < 0.05, na.rm = T)
sum(edgeR_results$table$FDR < 0.05, na.rm = T)
sum(nbsr_res$padj < 0.05, na.rm = T)

cols <- RColorBrewer::brewer.pal(8, "Set1")
RColorBrewer::brewer.pal(8, "Set1")
RColorBrewer::display.brewer.pal(8, "Set1")
cols2 <- c("#00AFBB", "#FC4E07", "#E7B800")
cols3 <- c("#56B4E9", "#E69F00")
cols5 <- c("black", "red")
#cols_methods <- c("Both"="#56B4E9", "NBSR"="#009E73", "EdgeR"="#E69F00", "DESeq2"="#CC79A7")

UPPER_THRESHOLD <- 0.05
NS_THRESHOLD <- 0.1

FEATURES_TO_PLOT <- 8
Y_OFFSET <- -0.006
LOG2FC_THRESHOLD <- 1

# Normalized counts showing the differing spread across the conditions.
deseq2_normalized_cts <- counts(dds2, normalized=TRUE)
X_dt <- as.data.table(X[,c("sample_id", "cell_tissue")])
X_dt$short_name <- gsub(pattern = "_lymphocyte_", " ", X_dt$cell_tissue)

# Get normalized counts.
cts_dt <- data.table(deseq2_normalized_cts,
                     features=rownames(Y))
cts_dt_melt <- melt(cts_dt, id.vars = "features", variable.name = "sample_id", value.name = "value")
cts_dt_melt_join <- dplyr::left_join(cts_dt_melt, X_dt, by = "sample_id")

# Props
props <- apply(Y, 2, normalize)
#props_dt <- data.table(t(props), sample_id=X$sample_id, cell_tissue=X$cell_tissue)
props_dt <- data.table(t(props), sample_id=X$sample_id)
props_dt_melt <- melt(props_dt, id.vars = c("sample_id"), variable.name = "features")
props_dt_melt_join <- dplyr::left_join(props_dt_melt, X_dt, by = "sample_id")

## 1. Make plots: normalized cts, sample props, and dispersion estimation.

# Uses objects in the global environment.
generate_dispersion_plots <- function(f_idxs, file_prefix="", file_suffix="", props=FALSE, plot_dispersion=TRUE)
{
  nbsr_pi1 <- nbsr_pi[f_idxs,celltype1_idx[1]]
  nbsr_pi2 <- nbsr_pi[f_idxs,celltype2_idx[1]]
  
  pi_dt1 <- data.table(feature = names(nbsr_pi1), pi=nbsr_pi1, short_name=X_dt$short_name[celltype1_idx[1]])
  pi_dt2 <- data.table(feature = names(nbsr_pi2), pi=nbsr_pi2, short_name=X_dt$short_name[celltype2_idx[1]])
  pi_dt <- rbind(pi_dt1, pi_dt2)

  nbsr_phi1 <- nbsr_phi[f_idxs,celltype1_idx]
  nbsr_phi2 <- nbsr_phi[f_idxs,celltype2_idx]

  phi_dt1 <- data.table(feature = rownames(Y)[f_idxs], 
                        edgeR_phi = edgeR_phi[f_idxs],
                        deseq2_phi = deseq2_phi[f_idxs],
                        nbsr_phi = nbsr_phi1)
  phi_dt1 <- melt(phi_dt1, id.vars = c("feature", "deseq2_phi", "edgeR_phi"), variable.name = "sample", value.name = "nbsr_phi")
  phi_dt2 <- data.table(feature = rownames(Y)[f_idxs], 
                       deseq2_phi = deseq2_phi[f_idxs],
                       edgeR_phi = edgeR_phi[f_idxs],
                       nbsr_phi = nbsr_phi2)
  phi_dt2 <- melt(phi_dt2, id.vars = c("feature", "deseq2_phi", "edgeR_phi"), variable.name = "sample", value.name = "nbsr_phi")
  phi_dt1$cell_type <- celltype1_short
  phi_dt2$cell_type <- celltype2_short
  phi_dt_plt <- rbind(phi_dt1, phi_dt2)
  phi_dt_plt[,feature_short := gsub("hsa-", "", feature)]
  phi_dt_plt$cell_type <- factor(phi_dt_plt$cell_type, levels = c(celltype1_short, celltype2_short))

  # Plot the dispersion across two conditions.
  pl <- ggplot(phi_dt_plt, aes(cell_type, nbsr_phi)) + theme_cowplot() +
    geom_boxplot(aes(fill="NBSR")) +
    geom_segment(aes(x = -Inf, xend = Inf, y = deseq2_phi, yend = deseq2_phi, color = 'DESeq2'), 
                 linetype = "dashed", linewidth = 1) +  
    geom_segment(aes(x = -Inf, xend = Inf, y = edgeR_phi, yend = edgeR_phi, color = 'EdgeR'), 
                 linetype = "dotted", linewidth = 1) +  
    facet_grid(~ feature_short) +
    scale_fill_manual(name="", values = c("NBSR" = cols[4])) +
    scale_color_manual(name = "", 
                       values = c("DESeq2" = cols[1], "EdgeR" = cols[2]), 
                       labels = c("DESeq2" = "DESeq2", "EdgeR" = "EdgeR"))
  pl <- pl + scale_y_log10()
  pl <- decorate_figure(pl, ylab="Dispersion (Log10)", xlab="", xtext_rotate_angle = -90)
  filename <- paste0(file_prefix, "biological_cv", file_suffix, ".pdf")
  ggsave(filename = paste0(subfigure_path, "/", filename), pl, width=16)
  
  fnames <- rownames(Y)[f_idxs]
  
  for (fname in fnames)
  {
    norm_cts_mean <- 
    pl <- ggplot(cts_dt_melt_join[features == fname], aes(short_name, value)) + 
      geom_jitter(aes(col=short_name)) + theme_minimal() # + facet_grid(~ feature_short)
    pl <- pl + theme(legend.position = "none")
    pl <- pl + scale_color_manual(values=cols3[1:2])
    pl <- pl + geom_point(data = cts_dt_melt_join[features == fname,.(norm_cts_avg=mean(value)),by=.(short_name)],
                          aes(short_name, norm_cts_avg, shape="Average of normalized cts"), size=4) +
      scale_shape_manual(values=c("Average of normalized cts"=3)) +
      guides(color="none") + theme(legend.position = "top", legend.title = element_blank())
    pl0 <- decorate_figure(pl, ylab = "Normalized counts", xlab="", xtext_rotate_angle = -90, vjust=0.5)
    
    pl <- ggplot(props_dt_melt_join[features == fname], aes(short_name, value)) + 
      geom_jitter(aes(col=short_name)) + theme_minimal() # + facet_grid(~ feature_short)
    pl <- pl + scale_color_manual(values=cols3[1:2])
    pl <- pl + geom_point(data = pi_dt[feature == fname],
                          aes(short_name, pi, shape="NBSR estimate"), size=4) +
      scale_shape_manual(values=c("NBSR estimate"=4)) +
      guides(color="none") + theme(legend.position = "top", legend.title = element_blank())
    pl1 <- decorate_figure(pl, ylab = "Sample proportions", xlab="", xtext_rotate_angle = -90, vjust=0.5)
    
    pl2 <- ggplot(phi_dt_plt[feature == fname], aes(cell_type, nbsr_phi)) + theme_minimal() +
      geom_boxplot(aes(fill="NBSR")) +
      geom_segment(aes(x = -Inf, xend = Inf, y = deseq2_phi, yend = deseq2_phi, color = 'DESeq2'), 
                   linetype = "dashed", linewidth = 1) +  
      geom_segment(aes(x = -Inf, xend = Inf, y = edgeR_phi, yend = edgeR_phi, color = 'EdgeR'), 
                   linetype = "dotted", linewidth = 1) +  
      scale_fill_manual(name="", values = c("NBSR" = cols[4])) +
      scale_color_manual(name = "", 
                         values = c("DESeq2" = cols[1], "EdgeR" = cols[2]), 
                         labels = c("DESeq2" = "DESeq2", "EdgeR" = "EdgeR")) +
      scale_y_continuous(position = "left")
    pl2 <- decorate_figure(pl2, ylab="Dispersion", xlab="", xtext_rotate_angle = -90,
                           vjust=0.5)
    
    fname_short <- gsub("hsa-", "", fname)
    if (!plot_dispersion) {
      # Plot props and normalized counts side-by-side: ignores props flag.
      pl <- pl0 + pl1
      pl <- pl + plot_annotation(title = fname_short, theme = theme(plot.title=element_text(size=20, hjust=0.5)))
      filename <- paste0(file_prefix, gsub("/", "|", fname_short), "_props_norm_cts", file_suffix, ".pdf")
      print(filename)
      ggsave(filename = paste0(subfigure_path, "/", filename), pl, width=6)
      next
    }
    if (props) {
      pl <- pl1 + pl2
      pl <- pl + plot_annotation(title = fname_short, theme = theme(plot.title=element_text(size=20, hjust=0.5)))
      filename <- paste0(file_prefix, gsub("/", "|", fname_short), "_props_dispersion", file_suffix, ".pdf")
      print(filename)
      ggsave(filename = paste0(subfigure_path, "/", filename), pl, width=6)
    } else {
      pl <- pl0 + pl2
      pl <- pl + plot_annotation(title = fname_short, theme = theme(plot.title=element_text(size=20, hjust=0.5)))
      filename <- paste0(file_prefix, gsub("/", "|", fname_short), "_norm_cts_dispersion", file_suffix, ".pdf")
      print(filename)
      ggsave(filename = paste0(subfigure_path, "/", filename), pl, width=6)
    }
  }
}

# Compare dispersion.
M <- log1p(rowMeans(deseq2_normalized_cts))
sparsity1 <- rowMeans(Y[,celltype1_idx] == 0)
sparsity2 <- rowMeans(Y[,celltype2_idx] == 0)
sparsity <- rowMeans(Y == 0)
phi_dt1 <- data.table(feature = rownames(Y), 
                      edgeR_phi = edgeR_phi,
                      deseq2_phi = deseq2_phi,
                      M = M,
                      sparsity = sparsity1,
                      nbsr_phi[,celltype1_idx])
phi_dt1 <- melt(phi_dt1, id.vars = c("feature", "deseq2_phi", "edgeR_phi", "M", "sparsity"), variable.name = "sample", value.name = "nbsr_phi")
phi_dt2 <- data.table(feature = rownames(Y), 
                      deseq2_phi = deseq2_phi,
                      edgeR_phi = edgeR_phi,
                      M = M,
                      sparsity = sparsity2,
                      nbsr_phi[,celltype2_idx])
phi_dt2 <- melt(phi_dt2, id.vars = c("feature", "deseq2_phi", "edgeR_phi", "M", "sparsity"), variable.name = "sample", value.name = "nbsr_phi")
phi_dt1$cell_type <- celltype1_short
phi_dt2$cell_type <- celltype2_short
phi_dt <- rbind(phi_dt1, phi_dt2)
phi_dt[,feature_short := gsub("hsa-", "", feature)]
phi_dt$cell_type <- factor(phi_dt$cell_type, levels = c(celltype1_short, celltype2_short))

phi_dt_mRNA <- data.table(feature = rownames(Y), 
                      edgeR_phi = edgeR_phi,
                      deseq2_phi = deseq2_phi,
                      M = M,
                      sparsity = sparsity)
phi_dt_mRNA[,feature_short := gsub("hsa-", "", feature)]
pts <- phi_dt_mRNA[M > 3 & deseq2_phi > 5]
pl <- ggplot(phi_dt_mRNA, aes(M, deseq2_phi)) + geom_point() +
  theme_bw() +
  ggrepel::geom_label_repel(data = pts,
                            aes(M, deseq2_phi, label = feature_short)) +
  geom_point(data=pts, aes(M, deseq2_phi), col="red")
pl <- decorate_figure(pl, xlab_text = "Log1p(M)", ylab_text = "DESeq2 dispersion")
ggsave(filename = paste0(figure_path, "/mean_expr_DESeq2_dispersion.pdf"), pl)

pl2 <- ggplot(phi_dt_mRNA, aes(sparsity, deseq2_phi)) + geom_point() + theme_bw() +
  ggrepel::geom_label_repel(data = pts,
                            aes(sparsity, deseq2_phi, label = feature_short)) +
  geom_point(data = pts, aes(sparsity, deseq2_phi), col="red")
pl2 <- decorate_figure(pl2, xlab_text = "Sparsity", ylab_text = "DESeq2 dispersion")
ggsave(filename = paste0(figure_path, "/sparsity_DESeq2_dispersion.pdf"), pl2)

cts_dt_melt_join[,feature_short := gsub("hsa-", "", features)]
cts_dt_melt_join$short_name <- factor(cts_dt_melt_join$short_name, levels = c(celltype1_short, celltype2_short))
props_dt_melt_join$short_name <- factor(props_dt_melt_join$short_name, levels = c(celltype1_short, celltype2_short))

pl3 <- ggplot(cts_dt_melt_join[features %in% pts$feature], aes(short_name, value, col=short_name)) + 
  theme_bw() + geom_jitter() + facet_grid(~ feature_short) + scale_color_manual(values = cols3) +
  theme(legend.title = element_blank(), strip.text = element_text(size=18), legend.text = element_text(size=14))
pl3 <- decorate_figure(pl3, xlab_text = "", ylab_text = "Normalized counts")
ggsave(filename = paste0(figure_path, "/DESeq2_normalized_counts_jitter_plot.pdf"), pl3)

deseq2_mu_dt <- data.table(features=rownames(deseq2_mu), deseq2_mu)
deseq2_mu_dt_melt <- melt(deseq2_mu_dt, id.vars = "features", variable.name = "sample_id", value.name = "predicted")
temp_dt <- left_join(cts_dt_melt_join, deseq2_mu_dt_melt, by = c("sample_id", "features"))
pl4 <- ggplot(temp_dt[features %in% pts$feature], aes(value, predicted, col=short_name)) + geom_point() + 
  theme_bw() + facet_wrap(~ feature_short, scales = "free") + 
  scale_color_manual(values = cols3) + 
  theme(legend.title = element_blank(), strip.text = element_text(size=18), legend.text = element_text(size=14))
pl4 <- decorate_figure(pl4, xlab_text = "Normalized counts", ylab_text = "DESeq2 predicted")
ggsave(filename = paste0(figure_path, "/DESeq2_normalized_counts_vs_predicted.pdf"), pl4, height = 6, width = 9)

# Features with large positive fold change.
# miR-146b-5p demonstrates the point that we are trying to make -- both are dense but on different scale of counts.
# Therefore, there should be different dispersion values.
f_idxs_pos <- order(nbsr_log2_fc, decreasing = T)[1:FEATURES_TO_PLOT]
generate_dispersion_plots(f_idxs_pos, file_prefix="pos_fc_")

# We should also demonstrate how the dispersion estimation looks 
# like when there aren't too much of a fold change
f_idxs_null <- which(abs(nbsr_log2_fc) < 0.1)[1:FEATURES_TO_PLOT]
generate_dispersion_plots(f_idxs_null, file_prefix="no_fc_")

# Features with large negative fold change.
f_idxs_neg <- order(nbsr_log2_fc, decreasing = F)[1:FEATURES_TO_PLOT]
generate_dispersion_plots(f_idxs_neg, file_prefix="neg_fc_")

# 2. Composition plots.
# Get top 30 features.
pi_A <- nbsr_pi[,celltype1_idx[1]]
pi_B <- nbsr_pi[,celltype2_idx[1]]

# Convert to data frames and add a condition label
dt_A <- data.table(features = names(pi_A), probability = pi_A, condition = celltype1_short)
dt_B <- data.table(features = names(pi_B), probability = pi_B, condition = celltype2_short)

dt_top30_A <- dt_A %>% top_n(30, probability) %>% arrange(desc(probability))
dt_top30_B <- dt_B %>% top_n(30, probability) %>% arrange(desc(probability))
unique_features <- unique(c(dt_top30_A$feature, dt_top30_B$feature))
length(intersect(dt_top30_A$feature, dt_top30_B$feature))

dt_combined <- bind_rows(dt_A[features %in% unique_features], dt_B[features %in% unique_features])
dt_combined$miRNA_short <- gsub("hsa-", "", dt_combined$features)
# Order the miRNAs by probability in A.
temp_A<- dt_combined[condition == celltype1_short]
temp_A <- temp_A[order(probability, decreasing = T)]
dt_combined$miRNA_short <- factor(dt_combined$miRNA_short, levels = temp_A$miRNA_short)
dt_combined[,condition := factor(condition, levels = c(celltype1_short, celltype2_short))]

# Create a grouped bar chart
pl <- ggplot(dt_combined, aes(x = miRNA_short, y = probability)) +
  geom_bar(aes(fill = condition), stat = "identity", position = position_dodge(), width = 0.9) +
  labs(x = "", y = "Estimated proportion") +
  theme_minimal() +
  guides(fill = guide_legend(title=NULL))
#filename <- paste0(out_path, "/figures/top30_miRNAs.pdf")
#ggsave(filename = filename, pl, width=16, height = 4)

# What is the total proportion accounted for by the selected 35 miRNAs?
dt_combined[,sum(probability),by=.(condition)]

joined_dt <- left_join(dt_combined, nbsr_res)
joined_dt$x <- as.numeric(joined_dt$miRNA_short)
joined_dt$y <- as.numeric(joined_dt$probability)
top30_estimated_prop <- joined_dt[,.(x, y=max(y), log2FoldChange, padj, condition),by=.(features)]
top30_estimated_prop <- top30_estimated_prop[condition == celltype1_short]
top30_estimated_prop$xbegin <- top30_estimated_prop$x - 2/5
top30_estimated_prop$xend <- top30_estimated_prop$x + 2/5
#top30_estimated_prop$ybegin <- top30_estimated_prop$yend <- top30_estimated_prop$y
top30_estimated_prop$ybegin <- 0
top30_estimated_prop$yend <- max(top30_estimated_prop$y)
get_pval_category <- function(pval) {
  ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01, "**",
                ifelse(pval < 0.05, "*", "NS")))
}
top30_estimated_prop$pval_cat <- get_pval_category(top30_estimated_prop$padj)
top30_estimated_prop$pval_cat <- factor(top30_estimated_prop$pval_cat, levels = c("***", "**", "*", "NS"))
top30_estimated_prop$significance <- abs(top30_estimated_prop$log2FoldChange) > LOG2FC_THRESHOLD

pl2 <- pl + geom_label(data=top30_estimated_prop,
                      aes(x=(xbegin + xend)/2,
                          y=ybegin+Y_OFFSET,
                          label=sprintf("%.2f", log2FoldChange),
                          col=significance), 
                      size=3,
                      inherit.aes = FALSE)
pl2 <- pl2 + annotate(geom = "text", x=min(top30_estimated_prop$x)-1.12, y=Y_OFFSET, label="Log2 FC")
pl2 <- pl2 + scale_y_continuous(breaks = c(0.05, 0.1, 0.15), limits = c(Y_OFFSET-0.001, max(top30_estimated_prop$yend)*1.08))
pl2 <- pl2 + coord_cartesian(clip="off")
pl2 <- pl2 + geom_segment(data=top30_estimated_prop, aes(x=xbegin, xend=xend, y=yend*1.03, yend=yend*1.03), inherit.aes = FALSE)
pl2 <- pl2 + geom_text(data=top30_estimated_prop, aes(x=(xbegin+xend)/2, y=yend*1.05, label=pval_cat), 
                       col = "black",
                       inherit.aes = FALSE)
pl2 <- pl2 + scale_fill_manual(values = cols3)
pl2 <- pl2 + scale_color_manual(values = cols5, guide = FALSE)
pl2 <- pl2 + ggtitle("30 miRNAs with highest estimated proportion in either cell type")
pl2 <- pl2 + theme(plot.title = element_text(size = 20))  # Enlarge title size
pl2 <- decorate_figure(pl2, xlab_text = "", ylab_text = "Estimated Proportion", xtext_rotate_angle = -45)
filename <- paste0(figure_path, "/top30_miRNAs_log2fc.pdf")
ggsave(filename = filename, pl2, width=18, height = 6)

f1 <- "hsa-miR-150-5p"
nbsr_pi[f1,celltype1_idx[1]] - nbsr_pi[f1,celltype2_idx[1]]
f1 <- "hsa-miR-92a-3p"
nbsr_pi[f1,celltype1_idx[1]] - nbsr_pi[f1,celltype2_idx[1]]

f1 <- "hsa-miR-146b-5p"
nbsr_pi[f1,celltype1_idx[1]] - nbsr_pi[f1,celltype2_idx[1]]
f1 <- "hsa-miR-151b/151a-5p"
nbsr_pi[f1,celltype1_idx[1]] - nbsr_pi[f1,celltype2_idx[1]]

# 3. Volcano plot
# Generic volcano plot.
library(EnhancedVolcano)
pl_deseq2 <- EnhancedVolcano(res2, rownames(res2), "log2FoldChange", "padj", 
                             pCutoff=UPPER_THRESHOLD, FCcutoff = LOG2FC_THRESHOLD,
                             title="DESeq2", subtitle = "", 
                             caption=paste0("Significant count: ", sum(res2$padj < UPPER_THRESHOLD & abs(res2$log2FoldChange) > LOG2FC_THRESHOLD, na.rm = T)),
                             selectLab = rownames(res2)[head(order(res2$padj, decreasing = F), FEATURES_TO_PLOT)])
pl_edgeR <- EnhancedVolcano(edgeR_results$table, rownames(edgeR_results$table), "logFC", "FDR", 
                            pCutoff=UPPER_THRESHOLD, FCcutoff = LOG2FC_THRESHOLD, 
                            title="EdgeR", subtitle = "", 
                            caption=paste0("Significant count: ", sum(edgeR_results$table$FDR < UPPER_THRESHOLD & abs(edgeR_results$table$logFC) > LOG2FC_THRESHOLD, na.rm = T)),
                            selectLab = rownames(edgeR_results$table)[head(order(edgeR_results$table$FDR, decreasing = F), FEATURES_TO_PLOT)])
pl_nbsr <- EnhancedVolcano(nbsr_res, nbsr_res$features, "log2FoldChange", "padj", pCutoff=UPPER_THRESHOLD, FCcutoff = LOG2FC_THRESHOLD,
                           title="NBSR", subtitle = "", 
                           caption=paste0("Significant count: ", sum(nbsr_res$padj < UPPER_THRESHOLD & abs(nbsr_res$log2FoldChange) > LOG2FC_THRESHOLD, na.rm = T)),
                           selectLab = nbsr_res$features[head(order(nbsr_res$padj, decreasing = F), FEATURES_TO_PLOT)])

ggsave(filename = paste0(figure_path, "/volcano_DESeq2.pdf"), pl_deseq2, height = 8, width = 10)
ggsave(filename = paste0(figure_path, "/volcano_edgeR.pdf"), pl_edgeR, height = 8, width = 10)
ggsave(filename = paste0(figure_path, "/volcano_nbsr.pdf"), pl_nbsr, height = 8, width = 10)

# 4. MA plots.

# Retrieve log2fc and FDRs for all three methods.
log2fc_dt <- data.table(features=rownames(Y), NBSR = nbsr_log2_fc, DESeq2=res2$log2FoldChange, EdgeR=edgeR_results$table$logFC)
adj_pval_dt <- data.table(features=rownames(Y), NBSR_pval = nbsr_res$padj, DESeq2_pval=res2$padj, EdgeR_pval=edgeR_results$table$FDR)

# We are interested in understanding the differences between
# NBSR and mRNA-seq methods when it comes to log2 fold change
# and their significances.
# Identify mRNAs whose significance differ between mRNA-seq methods 
# and NBSR. 

# Find all miRNAs where the methods differ.
nbsr_adj_pval_dt <- adj_pval_dt[NBSR_pval <= UPPER_THRESHOLD & DESeq2_pval > NS_THRESHOLD & EdgeR_pval > NS_THRESHOLD,]
mRNA_adj_pval_dt <- adj_pval_dt[NBSR_pval > NS_THRESHOLD & DESeq2_pval <= UPPER_THRESHOLD & EdgeR_pval <= UPPER_THRESHOLD,]
nbsr_adj_pval_dt[,method := "NBSR+/mRNA-"]
mRNA_adj_pval_dt[,method := "NBSR-/mRNA+"]
sub_adj_pval_dt <- rbind(nbsr_adj_pval_dt, mRNA_adj_pval_dt)
sub_log2fc_dt <- log2fc_dt[features %in% sub_adj_pval_dt$features]
sub_log2fc_dt <- left_join(sub_log2fc_dt, sub_adj_pval_dt)
dim(sub_log2fc_dt)

abs_diff_nbsr_pi <- abs(nbsr_pi[sub_log2fc_dt$features,celltype1_idx[1]] - nbsr_pi[sub_log2fc_dt$features,celltype2_idx[1]])
fnames <- c("hsa-miR-148b-3p", "hsa-miR-532-5p", "hsa-miR-142-5p")
abs_diff_nbsr_pi[fnames]

sub_log2fc_dt[features %in% fnames]

# Also include miRNAs that would have been missed by thresholding on log2fc.
#miRNAs_to_label30 <- top30_estimated_prop[abs(log2FoldChange) < LOG2FC_THRESHOLD & padj < UPPER_THRESHOLD,features]
#miRNAs_to_label_dt <- log2fc_dt[features %in% miRNAs_to_label30]
#miRNAs_to_label_dt <- left_join(miRNAs_to_label_dt, adj_pval_dt)
#miRNAs_to_label_dt$method <- "NBSR+/Threshold-"
#dim(miRNAs_to_label_dt)
#sub_log2fc_dt <- rbind(miRNAs_to_label_dt, sub_log2fc_dt)

# MA-plots.
avg_normalized_cts <- cts_dt_melt_join[,.(M=mean(value)),by=.(features)]
avg_normalized_cts[,short_name := gsub("hsa-", "", features)]
avg_normalized_cts_log2fc <- left_join(avg_normalized_cts, log2fc_dt)
sub_avg_normalized_cts_log2fc_dt <- left_join(sub_log2fc_dt, avg_normalized_cts)

names(cols2) <- sort(unique(sub_log2fc_dt$method))

generate_ma_pl <- function(method_, show_label=TRUE, show_legend=FALSE, fc_threshold=0)
{
  sub_dat <- sub_avg_normalized_cts_log2fc_dt[method %in% method_]
  ma_pl <- ggplot(avg_normalized_cts_log2fc, aes(log1p(M), NBSR)) + geom_point(alpha=0.1)
  ma_pl <- ma_pl + theme_minimal()
  ma_pl <- ma_pl + geom_point(data = sub_dat, 
                              aes(log1p(M), NBSR, col=method))
  ma_pl <- ma_pl + scale_color_manual(values = cols2[method_])
  if (fc_threshold > 0) {
    ma_pl <- ma_pl + geom_hline(yintercept = fc_threshold, linetype="dotted", col='red')
    ma_pl <- ma_pl + geom_hline(yintercept = -fc_threshold, linetype="dotted", col='red')
  }
  
  if (show_legend) {
    ma_pl <- ma_pl + theme(legend.title = element_blank())
    ma_pl <- ma_pl + theme(legend.position = c(0.85, 0.1))
  } else {
    ma_pl <- ma_pl + theme(legend.position = "none")
  }
  ma_pl <- ma_pl + ggtitle(paste0(celltype1_short, " vs ", celltype2_short))
  if (length(method_) == 1) {
    ma_pl <- ma_pl + labs(subtitle = method_)
  }
  if (show_label) {
    ma_pl <- ma_pl + ggrepel::geom_label_repel(data = sub_dat, 
                                               aes(log1p(M), NBSR, label=short_name),
                                               nudge_y = ifelse(sub_dat$NBSR < 0, -3, 3))
  }
  ma_pl <- decorate_figure(ma_pl, xlab_text = "Log1p(M)", 
                           ylab_text = "NBSR Log2FC")
  return(ma_pl)
}
method_ <- "NBSR-/mRNA+"
ma_pl1 <- generate_ma_pl(method_, fc_threshold = 1)
ggsave(filename = paste0(figure_path, "/MA_plot_NBSR1.pdf"), ma_pl1, width = 6)
method_ <- "NBSR+/mRNA-"
ma_pl2 <- generate_ma_pl(method_, fc_threshold = 1)
ggsave(filename = paste0(figure_path, "/MA_plot_NBSR2.pdf"), ma_pl2, width = 6)
# method_ <- "NBSR+/Threshold-"
# ma_pl3 <- generate_ma_pl(method_, show_label = FALSE)
# ggsave(filename = paste0(figure_path, "/MA_plot_NBSR3.pdf"), ma_pl3, width = 6)

methods <- sort(unique(sub_log2fc_dt$method))
sub_dat <- sub_avg_normalized_cts_log2fc_dt[method %in% methods & log1p(M) > 5]
ma_pl <- generate_ma_pl(methods, show_label = FALSE, show_legend = TRUE, fc_threshold = 1)
ma_pl <- ma_pl + geom_label_repel(data = sub_dat, 
                                  aes(log1p(M), NBSR, label=short_name),
                                  nudge_y = ifelse(sub_dat$NBSR < 0, -3, 3))
ggsave(filename = paste0(figure_path, "/MA_plot_NBSR.pdf"), ma_pl)

# Generate normalized counts plot for NBSR+/mRNA- and NBSR-/mRNA+
features <- sub_log2fc_dt[method == "NBSR+/mRNA-", features]
f_idxs_nbsr_pos <- which(names(nbsr_log2_fc) %in% features)
generate_dispersion_plots(f_idxs_nbsr_pos, file_prefix ="nbsr_pos_mRNA_neg_")
generate_dispersion_plots(f_idxs_nbsr_pos, file_prefix ="nbsr_pos_mRNA_neg_", props = TRUE)
generate_dispersion_plots(f_idxs_nbsr_pos, file_prefix ="nbsr_pos_mRNA_neg_", plot_dispersion = FALSE)

features <- sub_log2fc_dt[method == "NBSR-/mRNA+", features]
f_idxs_mRNA_pos <- which(names(nbsr_log2_fc) %in% features)
generate_dispersion_plots(f_idxs_mRNA_pos, file_prefix ="nbsr_neg_mRNA_pos_")
generate_dispersion_plots(f_idxs_mRNA_pos, file_prefix ="nbsr_neg_mRNA_pos_", props = TRUE)
generate_dispersion_plots(f_idxs_mRNA_pos, file_prefix ="nbsr_neg_mRNA_pos_", plot_dispersion = FALSE)

# Generate normalized counts and dispersion plot for features of interest.
# Specific features to focus based on the MA-plot and the bar plot.
features_to_label <- c("hsa-miR-16-5p",
                       "hsa-miR-92a-3p",
                       "hsa-miR-142-5p",
                       "hsa-miR-146b-5p",
                       "hsa-miR-148b-3p",
                       "hsa-miR-150-5p",
                       "hsa-miR-151b/151a-5p",
                       "hsa-miR-532-5p")
# hsa-miR-148b-3p was called -'ve by NBSR but +'ve by mRNA methods and
# it has relatively large expression (log1p(M) > 5).
f1 <- "hsa-miR-148b-3p"
nbsr_pi[f1,celltype1_idx[1]] - nbsr_pi[f1,celltype2_idx[1]]
log2(nbsr_pi[f1,celltype1_idx[1]]) - log2(nbsr_pi[f1,celltype2_idx[1]])

# 5. Compare log2FC estimates.
# Highlight the miRNAs that are selected for visualization of the bar plot.
top30_features <- top30_estimated_prop$features
top30_log2fc_dt <- log2fc_dt[features %in% top30_features]

features_to_label_log2fc_dt <- log2fc_dt[features %in% features_to_label]
features_to_label_log2fc_dt$short_name <- gsub("hsa-", "", features_to_label_log2fc_dt$features)

pl <- ggplot(log2fc_dt, aes(NBSR, DESeq2)) + geom_point(alpha=0.1) + theme_minimal()
pl <- pl + ggtitle(paste0(celltype1_short, " vs ", celltype2_short))  +
  labs(subtitle = "Comparison of Log2 FC (top 30 miRNAs from either cell type in red)")
pl <- pl + geom_abline(slope = 1, intercept = 0)
pl <- pl + geom_point(data = top30_log2fc_dt,
                      aes(NBSR, DESeq2), col="red")
pl <- pl + geom_label_repel(data = features_to_label_log2fc_dt,
                           aes(NBSR, DESeq2, label=short_name),
                           nudge_y = ifelse(features_to_label_log2fc_dt$NBSR < 0, 4, -4))
pl <- decorate_figure(pl, xlab_text = "NBSR", ylab_text = "DESeq2")
pl
ggsave(filename = paste0(figure_path, "/log2FC_NBSR_DESeq2.pdf"), pl)


pl <- ggplot(log2fc_dt, aes(NBSR, EdgeR)) + geom_point(alpha=0.1) + theme_minimal()
pl <- pl + geom_abline(slope = 1, intercept = 0)
pl <- pl + ggtitle(paste0(celltype1_short, " vs ", celltype2_short))  +
  labs(subtitle = "Comparison of Log2 FC (top 30 miRNAs from either cell type in red)")
pl <- pl + geom_point(data = top30_log2fc_dt,
                      aes(NBSR, EdgeR), col="red")
pl <- pl + geom_label_repel(data = features_to_label_log2fc_dt,
                            aes(NBSR, EdgeR, label=short_name),
                            nudge_y = ifelse(features_to_label_log2fc_dt$NBSR < 0, 4, -4))
pl <- decorate_figure(pl, xlab_text = "NBSR", ylab_text = "EdgeR")
ggsave(filename = paste0(figure_path, "/log2FC_NBSR_EdgeR.pdf"), pl)

pl <- ggplot(log2fc_dt, aes(DESeq2, EdgeR)) + geom_point(alpha=0.1) + theme_minimal()
pl <- pl + geom_abline(slope = 1, intercept = 0)
pl <- pl + ggtitle(paste0(celltype1_short, " vs ", celltype2_short))  +
  labs(subtitle = "Comparison of Log2 FC (top 30 miRNAs from either cell type in red)")
pl <- pl + geom_point(data = top30_log2fc_dt,
                      aes(DESeq2, EdgeR), col="red")
pl <- pl + geom_label_repel(data = features_to_label_log2fc_dt,
                            aes(DESeq2, EdgeR, label=short_name),
                            nudge_y = ifelse(features_to_label_log2fc_dt$NBSR < 0, 4, -4))
pl <- decorate_figure(pl, xlab_text = "DESeq2", ylab_text = "EdgeR")
ggsave(filename = paste0(figure_path, "/log2FC_DESeq2_EdgeR.pdf"), pl)

f_to_label_idxs <- which(names(nbsr_log2_fc) %in% features_to_label)
generate_dispersion_plots(f_to_label_idxs, file_prefix ="feature_to_label_", props = FALSE)
generate_dispersion_plots(f_to_label_idxs, file_prefix ="feature_to_label_", props = TRUE)



###################### MA plot for DESeq2 and EdgeR #######################
ma_pl <- ggplot(avg_normalized_cts_log2fc, aes(log1p(M), DESeq2)) + geom_point(alpha=0.1)
ma_pl <- ma_pl + theme_minimal()
ma_pl <- ma_pl + geom_point(data = sub_avg_normalized_cts_log2fc_dt, 
                            aes(log1p(M), DESeq2, col=method))
ma_pl <- ma_pl + scale_color_manual(values = cols2)
ma_pl <- ma_pl + theme(legend.title = element_blank())
ma_pl <- ma_pl + theme(legend.position = c(0.85, 0.12))
ma_pl <- ma_pl + ggtitle(paste0("MA plot ", celltype1_short, " vs ", celltype2_short))
ma_pl <- decorate_figure(ma_pl, xlab_text = "Log1p(M)", ylab_text = "DESeq2 Log2FC")
ggsave(filename = paste0(figure_path, "/MA_plot_DESeq2.pdf"), ma_pl)

ma_pl <- ggplot(avg_normalized_cts_log2fc, aes(log1p(M), EdgeR)) + geom_point(alpha=0.1)
ma_pl <- ma_pl + theme_minimal()
ma_pl <- ma_pl + geom_point(data = sub_avg_normalized_cts_log2fc_dt, 
                            aes(log1p(M), EdgeR, col=method))
ma_pl <- ma_pl + scale_color_manual(values = cols2)
ma_pl <- ma_pl + theme(legend.title = element_blank())
ma_pl <- ma_pl + theme(legend.position = c(0.7, 0.1))
ma_pl <- ma_pl + ggtitle(paste0("MA plot ", celltype1_short, " vs ", celltype2_short))
ma_pl <- decorate_figure(ma_pl, xlab_text = "Log1p(M)", 
                         ylab_text = "EdgeR Log2FC")
ggsave(filename = paste0(figure_path, "/MA_plot_EdgeR.pdf"), ma_pl)

