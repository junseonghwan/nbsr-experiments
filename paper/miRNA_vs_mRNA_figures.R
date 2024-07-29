rm(list=ls())
# Generate Figures for the paper.
library(data.table)
library(microRNAome)
library(cowplot)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(stringr)
source("paper/functions.R")
source("paper/dirichlet_estimation_functions.R")

MIN_SAMPLES <- 5
FEATURE_COUNT <- 10
MIN_READS <- 10

figure_path <- "paper/figures/miRNA_mRNA/"
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = T)
}

#### Load and process data
data("microRNAome")
load("data/rseData_01_06_2020.RData")

### recount2
irm <- which(colData(rseData_01_06_2020)$project == "DRP001797" &
               duplicated(colData(rseData_01_06_2020)$sample))
rseData_01_06_2020 <- rseData_01_06_2020[, -irm]

## remove samples with average read length > 200 -- pairend end sequecing with 100bps each.
irm <- which(rseData_01_06_2020$avg_read_length > 200)
rseData_01_06_2020 <- rseData_01_06_2020[, -irm]

# Select cell types with >= 10 samples.
unique_celltypes <- names(which(table(rseData_01_06_2020$celltype) >= 10)) # 7 of them
ind <- which(rseData_01_06_2020$celltype %in% unique_celltypes)
rseData_01_06_2020 <- rseData_01_06_2020[,ind]

# Eliminate rows that are all zeros.
feature_rm <- (rowSums(assay(rseData_01_06_2020) > 0) == 0)
rseData_01_06_2020 <- rseData_01_06_2020[!feature_rm,]

### microRNAome
# Remove rows that are zero everywhere.
icell <- which(microRNAome$sample_category == "cell_type")
dat <- microRNAome[ ,icell]
dat$cell_tissue <- gsub(pattern = "_", replacement = " ", x = dat$cell_tissue)

## filter samples with less than 100,000 total counts
i_rm <- which(colSums(assay(dat)) < 1e5)
dat <- dat[ ,-i_rm]

unique_celltypes <- names(which(table(dat$cell_tissue) >= 40))
ind <- which(dat$cell_tissue %in% unique_celltypes)
dat <- dat[,ind]

##### Prepare recount2 data for plotting.
# Check how many features are expressed.
expressed_features <- colSums(assay(rseData_01_06_2020) >= MIN_READS)
library_sizes <- colSums(assay(rseData_01_06_2020))
summary(expressed_features)
c(sum(expressed_features >= 5000), sum(expressed_features >= 10000), sum(expressed_features >= 20000))
rseData_01_06_2020$celltype <- gsub("([a-z])([A-Z])", "\\1 \\2", rseData_01_06_2020$celltype)
table(rseData_01_06_2020$celltype)

# For each sample, normalize and compute cumulative proportions.
mRNA_sorted_props <- apply(assay(rseData_01_06_2020), 2, get_sorted_props2)
colSums(mRNA_sorted_props)
ent <- apply(mRNA_sorted_props, 2, entropy)
feature_count <- dim(mRNA_sorted_props)[1]
max_ent <- compute_max_entropy(feature_count)
normalized_entropy <- ent / max_ent

feature1 <- apply(mRNA_sorted_props, 2, max)
boxplot(feature1)
recount2.dt <- data.table(cell_type=rseData_01_06_2020$celltype, feature1=feature1, normalized_entropy=normalized_entropy)

### Prepare microRNAome data for plotting
cts <- assay(dat)
features_rm <- rowSums(cts > 0) == 0
miRNA_sorted_props <- apply(cts[!features_rm,], 2, get_sorted_props2)
mean(abs(colSums(miRNA_sorted_props) - 1))
ent <- apply(miRNA_sorted_props, 2, entropy)
# Maximal entropy for each sample.
compute_max_entropy <- function(feature_count) {
  return(entropy(rep(1/feature_count, feature_count)))
}
feature_count <- dim(miRNA_sorted_props)[1]
max_ent <- compute_max_entropy(feature_count)
normalized_entropy <- ent / max_ent

feature1 <- apply(miRNA_sorted_props, 2, max)
boxplot(feature1)
microRNAome.dt <- data.table(cell_type=dat$cell_tissue, feature1=feature1, normalized_entropy=normalized_entropy)
table(microRNAome.dt$cell_type)

## Print the tables of sample counts for each cell type.
library(xtable)
xt <- xtable(recount2.dt[,.N,by=.(cell_type)], label = "supp-tbl3")
names(xt) <- c("Cell type", "N")
print(xt, file = "paper/tables/tbl_recount2.tex", include.rownames=FALSE)

xt <- xtable(microRNAome.dt[,.N,by=.(cell_type)], label = "supp-tbl2")
names(xt) <- c("Cell type", "N")
print(xt, file = "paper/tables/tbl_microRNAome_celltypes.tex", include.rownames=FALSE)

### Cumulative plot

# Compute the number of features needed to account for 90% of the reads for each cell type.
mRNA_cutoffs90 <- apply(assay(rseData_01_06_2020), 2, find_cutoff, p = 0.9)
miRNA_cutoffs90 <- apply(cts, 2, find_cutoff, p = 0.9)

mRNA_dt <- data.table(feature_count=mRNA_cutoffs90, celltype=rseData_01_06_2020$celltype)
miRNA_dt <- data.table(feature_count=miRNA_cutoffs90, celltype=dat$cell_tissue)

miRNA_pl <- ggplot(miRNA_dt, aes(celltype, feature_count)) + geom_boxplot() + theme_bw()
miRNA_pl <- miRNA_pl + ylab("Number of miRNAs to which\n90% of reads are mapped") + 
  xlab("") +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1))

mRNA_pl <- ggplot(mRNA_dt, aes(celltype, feature_count)) + geom_boxplot() + theme_bw()
mRNA_pl <- mRNA_pl + 
  ylab("Number of mRNAs to which\n90% of reads are mapped") + 
  xlab("") +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste0(figure_path, "/Fig1A.png"), plot = miRNA_pl, device = "png", width = 6, height = 6)
ggsave(filename = paste0(figure_path, "/Fig1B.png"), plot = mRNA_pl, device = "png", width = 6, height = 6)

pls1 <- ggarrange(miRNA_pl, mRNA_pl, heights = 1, widths = 1, 
                  nrow = 1, ncol = 2, align = "h")
ggsave(filename = paste0(figure_path, "/Fig1AB.png"), plot = pls1, device = "png", width = 11, height = 6)

#### Entropy figure
names(recount2.dt)
names(microRNAome.dt)
#recount2.dt$cell_type <- gsub("([a-z])([A-Z])", "\\1 \\2", recount2.dt$cell_type)

miRNA_pl <- ggplot(microRNAome.dt, aes(y=cell_type, x=normalized_entropy)) + geom_boxplot() + theme_bw()
#pl4 <- decorate_figure(pl4, xlab_text = "Cell Type", ylab_text = "Normalized entropy")
#pl4 <- decorate_figure(pl4, xlab_text = "", ylab_text = "Normalized entropy")
miRNA_pl <- miRNA_pl + scale_x_continuous(limits = c(0, 1)) +
  ylab("") + theme(axis.ticks.y = element_blank()) + xlab("") + 
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14))
miRNA_pl <- miRNA_pl + geom_segment(aes(x = 0, xend = 1, y=0.5, yend = 0.5), arrow =arrow(length = unit(0.02, "npc")))
miRNA_pl <- miRNA_pl + geom_text(aes(x = 0.85, y=0.65, label = "Increasing diversity"))
mRNA_pl <- ggplot(recount2.dt, aes(y=cell_type, x=normalized_entropy)) + geom_boxplot() + theme_bw()
#pl3 <- decorate_figure(pl3, xlab_text = "Cell Type", ylab_text = "Normalized entropy")
#pl3 <- decorate_figure(pl3, xlab_text = "", ylab_text = "Normalized entropy")
mRNA_pl <- mRNA_pl + scale_x_continuous(limits = c(0, 1)) + 
  theme(axis.ticks.y = element_blank()) + ylab("") + xlab("Normalized Entropy") + 
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14))
        
pls2 <- ggarrange(miRNA_pl, mRNA_pl, 
                   ncol = 1, nrow = 2, align = "hv")
ggsave(filename = paste0(figure_path, "/Fig1CD.png"), pls2, height = 5, width = 10)

### Tile figure

# Compute the IQR for each of the feature in the top features.
mRNA_dt <- data.table(cell_type = rseData_01_06_2020$celltype, t(mRNA_sorted_props[1:FEATURE_COUNT,]))
mRNA_dt <- melt(mRNA_dt, id.vars = "cell_type")
compute_iqr(mRNA_dt[cell_type == unique_celltypes[1] & variable == "V2",value])
mRNA_iqr_dt <- mRNA_dt[,.(IQR=compute_iqr(value), Median=median(value)), by=.(cell_type, variable)]
names(mRNA_iqr_dt) <- c("cell_type", "feature", "IQR", "Median")
mRNA_iqr_dt$data <- "mRNA"

miRNA_dt <- data.table(cell_type = dat$cell_tissue, t(miRNA_sorted_props[1:FEATURE_COUNT,]))
miRNA_dt <- melt(miRNA_dt, id.vars = "cell_type")
compute_iqr(miRNA_dt[cell_type == unique_celltypes[1] & variable == "V2",value])
miRNA_iqr_dt <- miRNA_dt[,.(IQR=compute_iqr(value), Median=median(value)), by=.(cell_type, variable)]
names(miRNA_iqr_dt) <- c("cell_type", "feature", "IQR", "Median")
miRNA_iqr_dt$data <- "miRNA"

tile_plot <- function(iqr_dt, limit_max)
{
  iqr_dt$feature <- as.numeric(str_replace(iqr_dt$feature, "V", ""))
  # Sort by feature 1.
  temp_dt <- iqr_dt[feature == 1,]
  sorted_celltypes <- temp_dt$cell_type[order(temp_dt$Median, decreasing = F)]
  iqr_dt$cell_type <- factor(iqr_dt$cell_type, levels = sorted_celltypes)
  pl <- ggplot(iqr_dt, aes(feature, cell_type)) + theme_bw()
  pl <- pl + geom_tile(aes(fill=Median, width = width, height = height))
  pl <- pl + xlab("Feature") + ylab("")
  pl <- pl + theme(axis.ticks = element_blank()) + scale_x_discrete(limits = as.character(1:FEATURE_COUNT))
  pl <- pl + scale_fill_viridis_c(limits = c(0, limit_max))
  #scale_fill_gradient(limits = c(0, 0.5), low = "blue", high = "red")
  return(pl)
}
limit_max <- max(c(miRNA_iqr_dt$Median, mRNA_iqr_dt$Median))
miRNA_iqr_dt$height <- miRNA_iqr_dt$width <- sqrt(miRNA_iqr_dt$Median) / sqrt(limit_max)
mRNA_iqr_dt$height <- mRNA_iqr_dt$width <- sqrt(mRNA_iqr_dt$Median) / sqrt(limit_max)

#mRNA_iqr_dt$cell_type <- gsub("([a-z])([A-Z])", "\\1 \\2", mRNA_iqr_dt$cell_type)
miRNA_pl <- tile_plot(miRNA_iqr_dt, limit_max) +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14))

mRNA_pl <- tile_plot(mRNA_iqr_dt, limit_max) +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14))

pls3 <- ggarrange(miRNA_pl, mRNA_pl, 
                ncol = 2, nrow = 1,
                common.legend = TRUE, legend = "bottom")
ggsave(filename = paste0(figure_path, "/Fig1EF.png"), pls3, height = 8, width = 16)

# The correlation figure (Figure 3 from the R01 proposal).
## proportion of non-zero counts
p_non_zero <- colSums(assay(dat)>0)

## proportion of reads going to each miRNA
p_reads <- sweep(assay(dat), 2, colSums(assay(dat)), FUN = "/")

mir_fig3 <- function(celltype, title){
  ind <- which(dat$cell_tissue == celltype)
  itop <- order(rowMeans(p_reads[ ,ind]), decreasing = TRUE)[1:2]
  df <- data.frame(miR1 = p_reads[itop[1], ind],
                   miR2 = p_reads[itop[2], ind])
  xycor <- cor(x=p_reads[itop[1], ind], 
               y=p_reads[itop[2], ind], 
               method="spearman")
  ggplot(df, aes(y=miR2, x=miR1)) + geom_point() +
    labs(tag = title, title = gsub("_", " ", celltype),
         subtitle = paste("Spearman Correlation =", round(xycor, 2))) +
    xlab(paste("Proportion of Reads Mapped to", 
               gsub("hsa-", "", rownames(p_reads)[itop[1]]))) + 
    ylab(paste("Proportion of Reads Mapped to", 
               gsub("hsa-", "", rownames(p_reads)[itop[2]]))) +  
    theme(plot.tag = element_text(size = 24, margin = margin(b = -20), 
                                  face = "bold")) +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
    theme(plot.subtitle = element_text(size = 14, hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 16)) + 
    theme(axis.title.y = element_text(size = 16)) +
    theme(axis.text = element_text(size=14)) +
    theme_cowplot()
}

fig3a <- mir_fig3("Endothelial", "")
fig3b <- mir_fig3("T lymphocyte CD8", "")
fig3c <- mir_fig3("Monocyte", "")

pls4 <- ggarrange(fig3a, fig3b, fig3c,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
ggsave(filename = paste0(figure_path, "/Fig1_correlation_miRNA.png"), pls4, height = 8, width = 16)

##### for recount data
## proportion of non-zero counts
p_non_zero <- colMeans(assay(rseData_01_06_2020) > 0)

## proportion of reads going to each mRNA
p_reads <- sweep(assay(rseData_01_06_2020), 2, 
                 colSums(assay(rseData_01_06_2020)), FUN = "/")

mrna_fig3 <- function(celltype, title){
  df <- data.frame(mR1 = p_reads[itop[1], ind],
                   mR2 = p_reads[itop[2], ind])
  xycor <- cor(x=p_reads[itop[1], ind], 
               y=p_reads[itop[2], ind], 
               method="spearman")
  ggplot(df, aes(y=mR2, x=mR1)) + geom_point() +
    labs(tag = title, title = gsub("_", " ", celltype),
         subtitle = paste("Spearman Correlation =", round(xycor, 2))) +
    xlab(paste("Proportion of Reads Mapped to", 
               rowData(rseData_01_06_2020)$symbol[itop[1]])) + 
    ylab(paste("Proportion of Reads Mapped to", 
               rowData(rseData_01_06_2020)$symbol[itop[2]])) + 
    theme(plot.tag = element_text(size = 24, margin = margin(b = -20), 
                                  face = "bold")) +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
    theme(plot.subtitle = element_text(size = 14, hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 16)) + 
    theme(axis.title.y = element_text(size = 16)) +
    theme(axis.text = element_text(size=14)) +
    theme_cowplot()
}

## Endothelial
ind <- which(colData(rseData_01_06_2020)$celltype == "Endothelial Cell")

## remove technical replicates from DRP001797
#irm <- which(colData(rseData_01_06_2020)$project[ind] == "DRP001797" &
#               duplicated(colData(rseData_01_06_2020)$sample[ind]))
#ind <- ind[-irm]

itop <- order(rowMeans(p_reads[ ,ind]), decreasing = TRUE)[1:20]
unlist(rowData(rseData_01_06_2020)$symbol[itop,])
## remove RNR2 (ribosome gene) and COX1 (mitochondrial gene)
itop <- itop[-c(1,3)]
fig3d <- mrna_fig3("EndothelialCell", "")

## Lymphocytes
ind <- which(colData(rseData_01_06_2020)$celltype == "Lymphocytes")
itop <- order(rowMeans(p_reads[ ,ind]), decreasing = TRUE)[1:20]
unlist(rowData(rseData_01_06_2020)$symbol[itop,])
## remove ribosome and mitochondrial genes and lncRNAs
itop <- itop[c(3,6,7)]
fig3e <- mrna_fig3("Lymphocytes", "")

## Macrophage
ind <- which(colData(rseData_01_06_2020)$celltype == "Macrophage")
itop <- order(rowMeans(p_reads[ ,ind]), decreasing = TRUE)[1:20]
unlist(rowData(rseData_01_06_2020)$symbol[itop,])
## remove ribosome and mitochondrial genes and lncRNAs
itop <- itop[c(1,3:5)]
fig3f <- mrna_fig3("Macrophage", "")

pls5 <- ggarrange(fig3d, fig3e, fig3f,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
ggsave(filename = paste0(figure_path, "/Fig1_correlation_mRNA.png"), pls5, height = 8, width = 16)

## For figure panel, we will use Endothelial cells.
pls6 <- ggarrange(fig3a, fig3b, 
                  ncol = 2, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
ggsave(filename = paste0(figure_path, "/Fig1GH.png"), pls6, height = 5, width = 10)

pls7 <- ggarrange(fig3d, fig3e, 
                  ncol = 2, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
ggsave(filename = paste0(figure_path, "/Fig1IJ.png"), pls7, height = 5, width = 10)

############################################################
# Supplementary figures
############################################################
# Compute the overall proportion taken up by top 10 features.
K <- 10
dim(miRNA_sorted_props[1:K,])
miRNA_top_expressed <- data.table(SampleName=dat$sample_id, CellType=dat$cell_tissue, t(miRNA_sorted_props[1:K,]))
mRNA_top_expressed <- data.table(SampleName=rseData_01_06_2020$sample, CellType=rseData_01_06_2020$celltype, t(mRNA_sorted_props[1:K,]))
names(miRNA_top_expressed) <- c(names(miRNA_top_expressed)[1:2], 1:K)
names(mRNA_top_expressed) <- c(names(mRNA_top_expressed)[1:2], 1:K)
miRNA_top_expressed <- melt(miRNA_top_expressed, id.vars = c("SampleName", "CellType"))
mRNA_top_expressed <- melt(mRNA_top_expressed, id.vars = c("SampleName", "CellType"))

pl <- ggplot(miRNA_top_expressed, aes(x=variable, y = value)) + geom_boxplot() + theme_bw()
pl <- pl + xlab("Feature") + ylab("Proportion expressed")
ggsave(filename = paste0(figure_path, "/miRNA_features", K, ".pdf"), pl)

pl <- ggplot(mRNA_top_expressed, aes(x=variable, y = value)) + geom_boxplot() + theme_bw()
pl <- pl + xlab("Feature") + ylab("Proportion expressed")
ggsave(filename = paste0(figure_path, "/mRNA_features", K, ".pdf"), pl)

# Feature 1 proportion
pl <- ggplot(microRNAome.dt, aes(x=cell_type, y=feature1)) + geom_boxplot() + theme_bw()
pl <- decorate_figure(pl, xlab_text = "Cell Type", ylab_text = "Proportion of reads taken by the most highly expressed feature")
pl <- pl + geom_hline(yintercept = mean(microRNAome.dt$feature1), col = 'red') + ggtitle(paste0("Overall mean: ", round(mean(microRNAome.dt$feature1), 3)))
ggsave(filename = paste0(figure_path, "/miRNA_feature1.pdf"), pl)

pl <- ggplot(recount2.dt, aes(x=cell_type, y=feature1)) + geom_boxplot() + theme_bw()
pl <- decorate_figure(pl, xlab_text = "Cell Type", ylab_text = "Proportion of reads taken by the most highly expressed feature")
pl <- pl + geom_hline(yintercept = mean(recount2.dt$feature1), col = 'red') + ggtitle(paste0("Overall mean: ", round(mean(recount2.dt$feature1), 3)))
ggsave(filename = paste0(figure_path, "/mRNA_feature1.pdf"), pl)

# Plot for percent miRNAs.
cdata <- as.data.table(colData(dat))
pl <- ggplot(cdata, aes(percent_miRNAs, log(total_input_reads))) + geom_point()
pl <- pl + theme_bw()
pl <- pl + stat_cor(method = "pearson", aes(label = paste(..r.label.., sep = "~`,`~")),
                    r.accuracy = 0.001, label.x = 0.01, label.y = 17)
pl <- pl + geom_smooth(method = "lm")
pl <- pl + xlab("% miRNAs") + ylab("Log of total input reads") + # + labs(tag="A") 
  theme(axis.title = element_text(size=20))
  
ggsave(filename = paste0(figure_path, "/percent_miRNA_vs_log_total_reads.pdf"), pl)

pl2 <- ggplot(cdata, aes(percent_miRNAs, log(filtered_miRNA_reads))) + geom_point() + theme_bw()
pl2 <- pl2 + stat_cor(method = "pearson", aes(label = paste(..r.label.., sep = "~`,`~")),
                    r.accuracy = 0.001, label.x = 0.01, label.y = 17)
pl2 <- pl2 + geom_smooth(method = "lm")
pl2 <- pl2 + xlab("% miRNAs") + ylab("Log of miRNA reads") + # + labs(tag="B")
  theme(axis.title = element_text(size=20))
ggsave(filename = paste0(figure_path, "/percent_miRNA_vs_log_miRNA_reads.pdf"), pl2)

# supp1 <- ggarrange(pl, pl2, 
#                    ncol = 2, nrow = 1)
# ggsave(filename = "../NBSR/paper/figures/percent_miRNAs_reads.png", supp1)
# ggsave(filename = "../NBSR/paper/figures/percent_miRNAs_reads.pdf", supp1)

pl <- ggplot(cdata, aes(percent_miRNAs, unique_miRNAs)) + geom_point() + theme_bw() +
  labs(x="% miRNAs", y="# of unique miRNAs") + 
  theme(axis.title = element_text(size=20)) +
  geom_smooth()
ggsave(filename = paste0(figure_path, "/percent_miRNA_vs_miRNA_captured.pdf"), pl)

pl <- ggplot(cdata, aes(filtered_miRNA_reads, unique_miRNAs)) + geom_point() + theme_bw() +
  labs(x="miRNA reads", y="# of unique miRNAs") + 
  theme(axis.title = element_text(size=20)) +
  geom_smooth()
ggsave(filename = paste0(figure_path, "/miRNA_reads_vs_miRNA_captured.pdf"), pl)

pl <- ggplot(cdata, aes(total_input_reads, unique_miRNAs)) + geom_point() + theme_bw() +
  labs(x="Input reads", y="# of unique miRNAs") + 
  theme(axis.title = element_text(size=20)) +
  geom_smooth()
ggsave(filename = paste0(figure_path, "/input_reads_vs_miRNA_captured.pdf"), pl)

pl <- ggplot(cdata, aes(log(total_input_reads), log(filtered_miRNA_reads))) + geom_point() + theme_bw() +
  labs(x="Input reads (Log)", y="miRNA reads (Log)") + 
  theme(axis.title = element_text(size=20)) +
  geom_smooth()
ggsave(filename = paste0(figure_path, "/input_reads_vs_miRNA_reads.pdf"), pl)
