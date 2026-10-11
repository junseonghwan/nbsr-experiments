# Repository for reproducing results related to NBSR paper for miRNA-seq

The scripts for reproducing the figures for [Improved differential expression analysis of miRNA-seq data by modeling competition to be counted](https://www.biorxiv.org/content/10.1101/2024.05.07.592964v2) can be found under `paper/`. The accompanying datasets are available on [Zenodo](https://doi.org/10.5281/zenodo.13127864).

## Pre-processing

- estimate_microRNAome.R
- dirichlet_estimation_functions.R
- functions.R

First, install the microRNAome package available on Bioconductor [available here](https://www.bioconductor.org/packages/3.16/data/experiment/html/microRNAome.html) and download recount2 data from [Zenodo](https://doi.org/10.5281/zenodo.13127864).
Also download mature miRNA sequences for human from [here](https://mirgenedb.org/fasta/hsa?mat=1) as `data/hsa.fas`.

Run `estimate_microRNAome.R` to obtain MLE for each of the cell types considered in the paper. The output will be written in `data/`.

## Figures comparing miRNA vs mRNA

- miRNA_vs_mRNA_figures.R

## Simulated data results

- generate_simulation_data.R
- batch_process_simul_results.R
- simulation_figures.R
- simulation_supp_figures.R

The validation data file names in Zenodo follow convention of `swap[0-3].tar.gz`, with `0-3` indicating the number of samples `n = 3, 5, 10, 20`. 
There is one additional data set named `swap2_EB.tar.gz`; this dataset contains results from running NBSR with Empirical Bayes dispersion estimation on the same dataset contained in `swap2.tar.gz`.
Create a new directory `data/validation` and uncompress the downloaded validation dataset there. 
The downloaded data already contains the outputs from NBSR. The script used for generating the plots can be found in `simulation_figures.R` and `simulation_supp_figures.R`.

The simulated data was generated using the script `generate_simulation_data.R`. 

### Latent simulation

The latent extension is `paper/generate_data_total_rna_imbalance_latent.R`. It uses
the existing alpha estimates in `data/microRNAome_dir_params.csv`, the feature
filter in `data/hsa.fas`, and real counts in `data/microRNAome_cd4_counts.rds`.
These inputs are included in this repository. Install the R packages
`data.table`, `seqinr`, and `MASS`, then run from a clone:

```sh
Rscript paper/generate_data_total_rna_imbalance_latent.R
```

The script estimates covariance from log proportions in the 42 real CD4 samples,
then draws `eta ~ Normal(log(alpha), Sigma)`, sets `pi = softmax(eta)`, and draws
multinomial counts. There is no additional Dirichlet draw or simulated z/gamma.
It uses the full covariance, with no PCA truncation or within-sample CLR step.
All-zero real features receive zero covariance. The pseudocount is `1e-6`, as
in the original alpha-estimation helper; it can strongly affect sparse features.

The LOW, MED, and HIGH replicates go to `data/total_imbalanced_logistic_normal/`,
under `<case>/sample_<n>/rep<r>/`. Each replicate contains `X.csv`, `Y.csv`,
`fc.csv`, and `eta.csv`; `Sigma.rds` is saved once at the output root, with miRNA
names on both dimensions. `Y.csv`, `fc.csv`, and `eta.csv` use the same retained
miRNAs. The settings remain 20 replicates and 10 samples per condition.
The old postprocessor that reads `z.csv` and `gamma.csv` needs adaptation to
use `eta.csv` for this generator.

See [CD4 data provenance](data/microRNAome_cd4_README.md) for the source and filters.
The bundled `hsa.fas` is the mature human miRNA FASTA from miRGeneDB linked
above.

## Immune cell analysis (T CD8 vs B CD19)

- immune_cells.R

Download the `immune.tar.gz` from Zenodo under `data/`. Then, decompress it, which should place data files and NBSR outputs under `data/immune/`.

## Colon adenocarcinoma cell line analysis

- carcinoma_generate_input.R
- carcinoma_analysis.R
- carcinoma_supp.R

Download the `carcinoma.tar.gz` from Zenodo under `data/`. Then, decompress it, which should place data files under `data/carcinoma/`.

## Contact

Please raise an issue if any figure cannot be reproduced or is missing from the paper, or if you have questions about running the scripts.
