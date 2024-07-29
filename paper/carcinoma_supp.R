rm(list=ls())
library(cowplot)
library(data.table)
library(DESeq2)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(MCMCpack)
library(microRNAome)
library(patchwork)
library(seqinr)
library(stringr)
library(VennDiagram)
source("paper/functions.R")
set.seed(1)

cellline_levels <- c("DLD-1", "DKO-1", "DKS-8")
data_path <- "data/carcinoma/"
figure_path <- paste0("paper/figures/carcinoma/")
if (!dir.exists(figure_path)) {
  dir.create(figure_path)
}

X <- read.csv(file = paste0(data_path, "/X.csv"))
Y <- as.matrix(read.csv(file = paste0(data_path, "/Y.csv")))
X$cell_tissue <- factor(X$cell_tissue, levels = cellline_levels)
dim(X)
dim(Y)

# For each cell line, generate plots as in Figure 1.
sorted_props <- apply(Y, 2, get_sorted_props2)
head(sorted_props)
ent <- apply(sorted_props, 2, entropy)
# Maximal entropy for each sample.
compute_max_entropy <- function(feature_count) {
  return(entropy(rep(1/feature_count, feature_count)))
}
feature_count <- dim(sorted_props)[1]
max_ent <- compute_max_entropy(feature_count)
normalized_entropy <- ent / max_ent

norm_ent_dt <- data.table(cell_type=X$cell_tissue, sample=X$sample_id, normalized_entropy=normalized_entropy)
miRNA_pl <- ggplot(norm_ent_dt, aes(y=cell_type, x=normalized_entropy)) + geom_boxplot() + theme_classic()
miRNA_pl <- miRNA_pl + scale_x_continuous(limits = c(0, 1)) + 
  ylab("") + theme(axis.ticks.y = element_blank()) + xlab("Shannon diversity (normalized)") + 
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14))
#miRNA_pl <- miRNA_pl+ coord_cartesian(ylim=c(0, 4))
miRNA_pl <- miRNA_pl + geom_segment(aes(x = 0, xend = 1, y=0.5, yend = 0.5), arrow =arrow(length = unit(0.02, "npc")))
miRNA_pl <- miRNA_pl + geom_text(aes(x = 0.85, y=0.65, label = "Increasing diversity"))
ggsave(filename = paste0(figure_path, "/shannon_diversity.pdf"), miRNA_pl)

cutoffs <- apply(Y, 2, find_cutoff, p = 0.90)
miRNA_dt <- data.table(feature_count=cutoffs, celltype=X$cell_tissue)
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
ggsave(filename = paste0(figure_path, "/feature_count90.pdf"), plot = miRNA_pl, width = 11, height = 6)

FEATURE_COUNT <- 15
miRNA_dt <- data.table(cell_type = X$cell_tissue, t(sorted_props[1:FEATURE_COUNT,]))
miRNA_dt <- melt(miRNA_dt, id.vars = "cell_type")
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
  pl <- pl + geom_tile(aes(fill=Median, width = width, height = height), col="black")
  #pl <- pl + geom_tile(aes(fill=Median), col="gray")
  pl <- pl + xlab("Feature") + ylab("")
  pl <- pl + theme(axis.ticks = element_blank()) + scale_x_discrete(limits = as.character(1:FEATURE_COUNT))
  pl <- pl + scale_fill_viridis_c(limits = c(0, limit_max))
  #scale_fill_gradient(limits = c(0, 0.5), low = "blue", high = "red")
  return(pl)
}
limit_max <- max(miRNA_iqr_dt$Median)
miRNA_iqr_dt$height <- miRNA_iqr_dt$width <- sqrt(miRNA_iqr_dt$Median) / sqrt(limit_max)

#mRNA_iqr_dt$cell_type <- gsub("([a-z])([A-Z])", "\\1 \\2", mRNA_iqr_dt$cell_type)
miRNA_pl <- tile_plot(miRNA_iqr_dt, limit_max) +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14))
ggsave(filename = paste0(figure_path, "/feature_prop_heatmap.pdf"), miRNA_pl, height = 2, width = 6)

# Bar plot of the median.
miRNA_iqr_dt$feature <- factor(miRNA_iqr_dt$feature, levels = rev(levels(miRNA_iqr_dt$feature)))
pl <- ggplot(miRNA_iqr_dt, aes(Median, cell_type, fill = Median)) + geom_bar(stat="identity", col="white")
pl <- pl + scale_fill_viridis_c(limits = c(0, limit_max))
pl <- pl + theme_cowplot()
pl <- pl + scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = seq(0, 0.8, 0.1))
pl <- pl + theme(axis.ticks.y = element_blank())
pl <- decorate_figure(pl, xlab_text = "Median expression", ylab_text = "")
ggsave(filename = paste0(figure_path, "/feature_prop_barplot.pdf"), pl, width=8, height=1.5)

cor(t(sorted_props[1:5,]), method = "spearman")

miRNA_iqr_dt[,sum(Median),by=.(cell_type)]
