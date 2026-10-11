# CD4 counts for covariance estimation

`microRNAome_cd4_counts.rds` is a plain numeric count matrix with 2,546 miRNA
rows and 42 sample columns. `microRNAome_cd4_samples.csv` contains matching
sample metadata; `count_column` matches the matrix column names.

Source: the `microRNAome` object (not `microRNAome_v2`) in microRNAome 1.20.0,
Bioconductor 3.16:
https://bioconductor.org/packages/3.16/data/experiment/src/contrib/microRNAome_1.20.0.tar.gz

The filters reproduce `paper/estimate_microRNAome.R` for CD4 cells:

- `sample_category == "cell_type"`
- `study_id == 89`
- `tolower(cell_tissue) == "t_lymphocyte_cd4"`
- at least 100,000 total counts, summed over the original 2,546 miRNAs

Counts have not been normalized or smoothed. All original miRNAs are retained
in this input file. The generator matches and filters them using the existing
alpha CSV and `hsa.fas`, leaving 596 miRNAs. Of these, 199 have zero counts in
all 42 real samples. Their covariance rows and columns are set to zero.

Re-estimating alpha on the full 2,546-by-42 matrix with the repository's
`estimate_dirichlet()` and 1,000 iterations reproduces the CD4 alpha CSV to
numerical precision (maximum absolute difference about 1.33e-9).
Newer microRNAome releases contain a different data object, so this source
version is fixed for reproducibility.

To extract these files again, download and unpack that archive, load
`microRNAome/data/microRNAome.rda` with the SummarizedExperiment package loaded,
apply the filters above, and save `assay(microRNAome, "counts")[, keep]` with
`saveRDS()`. The simulation itself does not require microRNAome or
SummarizedExperiment to be installed.
