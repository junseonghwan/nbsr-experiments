## grabbing some miRNA-seq from TCGA
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install('GenomicDataCommons')

library(GenomicDataCommons)
GenomicDataCommons::status()

# Step 1: Get case IDs for kidney samples
kidney_cases <- cases() |>
  filter(~ project.project_id == 'TCGA-KIRC' &
           diagnoses.tissue_or_organ_of_origin == 'kidney, nos' &
           demographic.race == "white" &
           demographic.ethnicity == "not hispanic or latino") |>
  results_all()

kidney_case_ids <- kidney_cases$case_id

# Step 2: Get files for those cases
kidney_files <- files() |>
  filter(~ cases.case_id %in% kidney_case_ids &
           data_type == 'miRNA Expression Quantification') |>
  results_all()

# Extract file IDs
file_ids <- kidney_files$file_id

# Check how many files
length(file_ids)

# Download to specific directory
fnames <- gdcdata(file_ids, destination_dir = "~/data/tcga_mirna_kidney")

# miRNA files are typically tab-delimited text files
# Read one file as an example
mirna_data <- read.table(fnames[1], header = TRUE, sep = "\t")
head(mirna_data)

# Or read all files into a list
all_mirna <- lapply(fnames, function(f) {
  read.table(f, header = TRUE, sep = "\t")
})

## test if all miRNA_IDs are identical
all(sapply(all_mirna, function(df) {
  all(df$miRNA_ID == all_mirna[[1]]$miRNA_ID)
}))

## merge counts
mirna_counts <- do.call(cbind, lapply(all_mirna, function(df) df$read_count))
rownames(mirna_counts) <- all_mirna[[1]]$miRNA_ID

## create SummarizedExperiment object
tcga_kirc_se = SummarizedExperiment(
  assays = list(counts = mirna_counts),
  rowData = DataFrame(miRNA_ID = all_mirna[[1]]$miRNA_ID),
  colData = DataFrame(file_id = colnames(mirna_counts))
)

# Query GDC for metadata
file_metadata <- files() |>
  filter(~ file_id %in% tcga_kirc_se$file_id) |>
  select(c('file_id',
           'cases.submitter_id',
           'cases.samples.sample_type',
           'cases.samples.submitter_id',
           'cases.diagnoses.ajcc_pathologic_stage')) |>
  results_all()

# Process nested structure
library(tidyr)
metadata_clean <- as_tibble(file_metadata) %>%
  unnest(cases, names_sep = ".") %>%
  unnest(cases.samples, names_sep = ".", keep_empty = TRUE) %>%
  unnest(cases.diagnoses, names_sep = ".", keep_empty = TRUE) %>%
  dplyr::select(file_id,
         case_id = cases.submitter_id,
         sample_type = cases.samples.sample_type,
         sample_id = cases.samples.submitter_id,
         ajcc_stage = cases.diagnoses.ajcc_pathologic_stage) %>%
  dplyr::distinct(file_id, .keep_all = TRUE)  # Remove duplicates if any

# Merge with existing colData
coldata_df <- as.data.frame(colData(tcga_kirc_se))
coldata_merged <- dplyr::left_join(coldata_df, metadata_clean, by = "file_id")

# Update colData in SummarizedExperiment
colData(tcga_kirc_se) <- DataFrame(coldata_merged)

save(tcga_kirc_se, file = "~/data/tcga_mirna_kidney/tcga_kirc_mirna_se.rda", compress = "xz")
