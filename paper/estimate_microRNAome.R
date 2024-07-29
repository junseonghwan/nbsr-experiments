rm(list=ls())
library(data.table)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(magrittr)
library(microRNAome)
library(MCMCpack)
library(tidyverse)
library(xtable)

source("paper/dirichlet_estimation_functions.R")
source("paper/functions.R")

seed <- 1
set.seed(seed)

# Keep cell_tissue with at least MIN_SAMPLES.
MIN_SAMPLES <- 20

data("microRNAome")
microRNAome$cell_tissue <- tolower(microRNAome$cell_tissue)
cell_types <- unique(microRNAome$cell_tissue)

cdat <- as.data.table(colData(microRNAome))
sample_categories <- unique(cdat$sample_category)
sample_counts <- cdat[,.N,by=.(sample_category,cell_tissue)]
# Most cancer cell lines have very few samples.
sample_counts[sample_category == "cancer_cell_line"]
sample_counts[sample_category == "tissue"]
sample_counts[N >= MIN_SAMPLES]

# Let's focus on cell_type.
icell <- which(microRNAome$sample_category == "cell_type")
dat <- microRNAome[ ,icell]
dat$cell_type <- gsub(pattern = "_", replacement = " ", x = dat$cell_tissue)

table(dat$cell_type, dat$study_id)

## filter samples with less than 100,000 total counts
i_rm <- which(colSums(assay(dat)) < 1e5)
dat <- dat[ ,-i_rm]

## We will focus on study id 89 to ensure uniformity in the samples.
table(dat$cell_tissue, dat$study_id)
table(dat$study_id)

ind <- which(dat$study_id == 89)
dat <- dat[,ind]
dat

#saveRDS(dat, file = "data/microRNAome_polya.RDS")

unique_celltypes <- unique(dat$cell_tissue)

table(dat$cell_tissue)
tbl <- data.frame(table(dat$cell_tissue))
tbl$Var1 <- stringr::str_to_title(tbl$Var1)
tbl$Var1 <- gsub("_", " ", tbl$Var1)
tbl$Var1 <- gsub("cd", "CD", tbl$Var1)
colnames(tbl) <- c("Cell type", "N")
tbl$CellType <- as.character(tbl$CellType)
#tbl$CellType <- tools::toTitleCase(tbl$CellType)
#tbl$CellType <- gsub(pattern = "_", replacement = " ", x = tbl$CellType)
sink(file = "paper/tables/tbl_microRNAome_polya.tex")
print(xtable(tbl, label = "supp-tbl1", 
             caption = "Sample counts for each cell type used in the microRNAome dataset analysis."),
      include.rownames = F)
sink()
sum(tbl$Count)

cts <- assays(dat)$counts
dim(cts)
# Keep miRNAs expressed in at least 3 samples.
#miRNA_list <- rownames(cts)[rowSums(cts > 0) >= 3]
miRNA_list <- rownames(cts)[rowSums(cts > 0) >= 0]
length(miRNA_list)

# Estimate Dirichlet parameters for each sample.
K <- length(miRNA_list)
print(paste0("Num miRNAs: ", K))
initial_alpha <- rep(1,K)
dirichlet_params <- data.frame()
sample_means <- list()
concentration_params <- list()
for (cell_tissue in unique_celltypes)
{
  print(cell_tissue)
  cell_tissue_col_idx <- which(grepl(paste0("^", cell_tissue, "$"), dat$cell_tissue))

  cell_type_dat <- dat[miRNA_list,cell_tissue_col_idx]
  cts <- assays(cell_type_dat)$counts

  sample_count <- dim(cts)[2]
  library_size <- colSums(cts)
  
  props <- apply(cts, 2, normalize)
  #sorted_props <- apply(cts, 2, get_sorted_props2)
  sample_means[[cell_tissue]] <- rowMeans(props)
  
  s_init <- 1
  ind <- which(sample_means[[cell_tissue]] > 0)
  s_new <- optimize_concentration(cts[ind,], sample_means[[cell_tissue]][ind], s_init, max_iter = 10000)
  concentration_params[[cell_tissue]] <- s_new$s
  #alpha_hat <- alpha_bar * s_new$s
  #head(sort(alpha_hat, decreasing = T), 10)

  ret <- estimate_dirichlet(t(cts), initial_alpha, max_iters = 1000, verbose = FALSE)
  print(ret$converged)
  alpha_hat <- ret$alpha_hat
  sum(alpha_hat)
  sum(alpha_hat > 1)
  sum(normalize(alpha_hat) > 0)
  
  dirichlet_params <- rbind(dirichlet_params, 
                            data.frame(cell_type = cell_tissue,
                                       alpha_hat = alpha_hat,
                                       miRNA = rownames(cts)))
}
write.table(dirichlet_params, file = "data/microRNAome_dir_params.csv", row.names = F, quote = F, sep=",")

dirichlet_params <- fread(file = "data/microRNAome_dir_params.csv")
conc_params <- dirichlet_params[,.(s=sum(alpha_hat)),by=.(cell_type)]

dirichlet_params[,.(s=sum(alpha_hat > 1)),by=.(cell_type)]
dirichlet_params[,.(s=sum(alpha_hat <= 1)),by=.(cell_type)]

