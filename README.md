# Repository for reproducing results related to NBSR paper for miRNA-seq

The scripts for reproducing the figures for [Improved differential expression analysis of miRNA-seq data by modeling competition to be counted](https://www.biorxiv.org/content/10.1101/2024.05.07.592964v2). The scripts for reproducing the figures can be found under `paper/`.

## Pre-processing

- estimate_microRNAome.R
- dirichlet_estimation_functions.R
- functions.R

As a first step, make sure to install microRNAome package available on [bioconductor](https://www.bioconductor.org/packages/3.16/data/experiment/html/microRNAome.html) and recount2 data from [Zenodo](https://doi.org/10.5281/zenodo.13127864).
Also download mature miRNA sequences for human from [here](https://mirgenedb.org/fasta/hsa?mat=1) and place it as `data/hsa.fas`.

Run `estimate_microRNAome.R` to obtain MLE for each of the cell types considered in the paper.

## Figures comparing miRNA vs mRNA

- miRNA_vs_mRNA_figures.R

## Simulated data results

- generate_simulation_data.R
- batch_process_simul_results.R
- simulation_figures.R
- simulation_supp_figures.R

First download the validation dataset from [Zenodo](https://doi.org/10.5281/zenodo.13127864). 
The file names follow convention of `swap[0-3].tar.gz`, with `0-3` indicating the number of samples `n = 3, 5, 10, 20`. 
There is one additional data set named `swap2_EB.tar.gz`; this dataset contains results from running NBSR with Empirical Bayes dispersion estimation on the same dataset contained in `swap2.tar.gz`.
Create directory `data/validation` and uncompress the downloaded validation dataset there. 
The downloaded data already contains the outputs from NBSR and the plots. The script used for generating the plots can be found in `simulation_figures.R` and `simulation_supp_figures.R`.

The simulated data was generated using script `generate_simulation_data.R`. 

## Immune cell analysis (T CD8 vs B CD19)

- immune_cells.R

## Colon adenocarcinoma cell line analysis

- carcinoma_generate_input.R
- carcinoma_analysis.R
- carcinoma_supp.R

## Contact

Please raise an issue if any figure cannot be reproduced/missing from the paper or if you have questions about running the scripts.

