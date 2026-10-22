# Reproducing "Bayesian compositional differential abundance with negative binomial softmax regression" (RECOMB 2026 submission)

Everything the paper reports comes from fits of NBSR (<https://github.com/junseonghwan/nbsr>) and of
six published methods to simulated and real sequencing data. This directory holds what is needed to
rebuild the paper's tables and figures from the released per-feature outputs, to regenerate the
simulated inputs, and to refit the models.

## 1. Rebuild the tables and figures from the released outputs

Download the outputs deposit from Zenodo (DOI to be added) and unpack it as `recomb2026/data/outputs/`.
It holds one gzipped CSV per benchmark with one row per method, simulation case, sample size,
replicate and feature, unfiltered (every feature a method reported, with a `tested` flag), plus a
README with the column dictionary and `SHA256SUMS`. Then, from the repository root:

```bash
Rscript recomb2026/scripts/paper_tables.R          # -> recomb2026/tables/*.csv and tables.md
Rscript recomb2026/scripts/make_paper_figures.R    # -> recomb2026/figures/
```

Both need R with `data.table`, `ggplot2` and `ragg`, and `gzip` on the path. `paper_tables.R`
writes Tables 1 to 3 and the appendix tables that are computed from per-feature summaries (the
calibration, decision-rule, Laplace-against-HMC and Wald-rule tables, the microbiome grid with power
at matched FDR, the per-sample dispersion variants, the permutation null, the real-label analyses and
the ablation variants). `make_paper_figures.R` writes the precision-recall and average-precision
panels of both benchmarks, the adjusted-p and lfsr scatter plots, the three calibration figures and
the concentration figure. Four items need the posterior draws themselves, which are not deposited
(about 100 MB per fit): the variance decomposition of the shift, the reference-set sensitivity
tables, the dispersion-regularisation factorial on the fitted composition, and Figure 2; the
scripts that produced them are in the private pipeline repository and the draws are available on
request.

## 2. Regenerate the simulated inputs

The inputs deposit (`paper_inputs.tar.gz`, DOI to be added) holds the exact replicates used, each as
`Y.csv` (counts), `X.csv` (design) and `fc.csv` (truth). They can also be regenerated:

- **miRNA benchmark**: the Dirichlet-multinomial generator of the miRNA paper,
  `../paper/generate_simulation_data.R` in this repository, seeded from the microRNAome CD4+
  T-lymphocyte estimate (`../data/`), at n = 10, 20, 40 per condition; cases LOW, MED and HIGH
  perturb 20, 10 and 5 features.
- **microbiome benchmark**: `generators/simulate_lndm.R`, the logistic-normal Dirichlet-multinomial
  generator anchored on the Baxter et al. (2016) healthy stool template in `data/templates/`
  (`estimate_dirichlet.R` and `estimate_lognormal.R` built the template files from MicrobiomeHD,
  which `prepare_microbiomehd.py` downloads); cases in `generators/simulation_cases.csv`, sample
  sizes 40 to 320. `--sigma 0` gives the plain Dirichlet-multinomial of the concentration sweep,
  `--A` sets the concentration, `--A_logsd` draws one per sample. `validate_generator.R` is the
  template comparison of the appendix.
- **SELEX**: `generators/prepare_selex.R` writes the ALDEx2 `selex` data (`data/selex/`) in the
  pipeline layout with the 27 characterised variants and the 20 random 5 + 5 subsets.
- **label permutations and real labels**: `prepare_permutation.R` and `prepare_real.R` on the
  processed Baxter and Schubert genus tables in `data/templates/`.

## 3. Refit the models

`model/` holds the Stan programs and the settings, so the fits can be repeated with any Stan
front end:

- `nbsr.stan`: the model of the paper. `methods.json` records every prior and sampler setting
  (four chains, 1,000 warm-up and 4,000 sampling iterations, `adapt_delta` 0.95, maximum tree
  depth 12; `sigma_beta2_intercept = 100` is the diffuse intercept prior used on the microbiome
  benchmark). Chains start at the prior centre with a small jitter, the dispersion intercept at a
  method-of-moments value.
- `nbsr_oracle.stan`, `nbsr_offset.stan`, `nbsr_disp.stan`, `nbsr_hs.stan`, `nbsr_nobj.stan` with
  `methods.ablation.json`: the appendix ablations.
- `estimand.R`: the per-draw estimand, the CLR effect re-centred on the mode of a Gaussian kernel
  density estimate over the reference set (the most abundant tenth of the panel with a floor of 30
  features on the miRNA benchmark; the whole panel elsewhere), and the summaries written per fit.
- `fit_laplace.py`: the Laplace route, NBSR's posterior mode and Gaussian draws from the Hessian,
  using the NBSR package (commit 67e11cb); `methods.selex.json` adds the SELEX prior scales
  elicited by its empirical-Bayes procedure.

The competitors were run at their package defaults as the paper's appendix states
(ANCOMBC 2.14.0, ALDEx2 1.42.0, MicrobiomeStat 1.4, DESeq2 1.50.2, edgeR 4.8.2, apeglm 1.32.0);
the exact calls are in the appendix "Competitor methods".

The Slurm pipeline that ran all of this on the Empire AI cluster is specific to that system and is
not released; `scripts/export_paper_data.R` is the script that produced the deposit from its
per-fit summaries.
