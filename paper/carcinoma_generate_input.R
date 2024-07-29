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
cols4 <- c("red", "blue")
cols5 <- c("black", "red")
cols_cell_line <- c("DKO-1"="#009E73", "DKS-8"="#CC79A7", "DLD-1"="#F0E442")

out_path <- "data/carcinoma/"
WRITE_DESeq2_DISPERSION <- FALSE
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = T)
}

data("microRNAome")

### microRNAome
sort(table(microRNAome[,microRNAome$sample_category == "cancer_cell_line"]$cell_tissue))
sort(table(microRNAome$organ))

carcinoma_type <- "Colon adenocarcinoma"
sample_i <- grepl(pattern = carcinoma_type, microRNAome$organ)
sum(sample_i)
dat <- microRNAome[,sample_i]
table(dat$sample_id)
i_rm <- (colSums(assay(dat)) < 1e6)
dat <- dat[,!i_rm]
dim(dat)
table(dat$study_id, dat$cell_tissue)
table(dat$study_id, dat$organ)
table(dat$organ, dat$cell_tissue)
i_keep <- (dat$study_id %in% c(17))
dat <- dat[,i_keep]
dim(dat)

cts <- assay(dat)

#Y <- cts[!(rowMeans(cts == 0) == 1),]
MIN_SAMPLES <- 3
Y <- cts[(rowSums(cts > 0) >= MIN_SAMPLES),]
X <- colData(dat)

X$cell_tissue <- toupper(gsub("([a-z]+)([0-9]+)", "\\1-\\2", X$cell_tissue))
celltype1 <- "DLD-1"
celltype2 <- "DKO-1"
celltype3 <- "DKS-8"
cellline_levels <- c(celltype1, celltype2, celltype3)
X$cell_tissue <- factor(X$cell_tissue, levels = cellline_levels)
dim(Y)
dim(X)

# We keep more than 99.9% of the reads with the above filtering.
colSums(Y)/colSums(cts)

write.table(X, file = paste0(out_path, "/X.csv"), quote = F, col.names = T, row.names = F, sep=",")
write.table(Y, file = paste0(out_path, "/Y.csv"), quote = F, col.names = T, row.names = T, sep=",")

# Run DESeq2 and output deseq2_mu.csv
se <- SummarizedExperiment(assays = list(counts = Y), colData = X)
dds2 <- DESeqDataSet(se, ~ cell_tissue)
dds2 <- DESeq(dds2, fitType = "local")
rdata <- rowData(dds2)
deseq2_mu <- assays(dds2)[["mu"]]
deseq2_phi <- rdata$dispersion
write.table(deseq2_mu, file = paste0(out_path, "/deseq2_mu.csv"), quote=F, col.names = T, row.names=F, sep=",")
if (WRITE_DESeq2_DISPERSION) {
  names(deseq2_phi) <- rownames(Y)
  write.table(deseq2_phi, file = paste0(out_path, "/dispersion.csv"), quote=F, col.names = F, row.names=F, sep=",")
}
