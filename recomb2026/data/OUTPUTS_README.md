# Per-feature outputs of every fit in the paper

One row per method, simulation case, sample size (per condition), replicate and
feature, **unfiltered**: every feature a method's summary carries, no restriction
to a common set (reapply one with `tested`). Produced by
`scripts/export_paper_data.R` in nbsr-benchmarks from each unit's `summary.rds`;
`scripts/make_paper_figures.R` here rebuilds the paper's benchmark figures from
these files alone (`Rscript scripts/make_paper_figures.R` from the paper
directory; needs data.table, ggplot2, ragg and `gzip` on the path).

| file | benchmark | methods |
|---|---|---|
| `mirna_nbsr_features.csv.gz`, `mirna_compositional_features.csv.gz`, `mirna_countbased_features.csv.gz` | miRNA (LOW/MED/HIGH; n = 10, 20, 40; 20 replicates) | NBSR by HMC and Laplace at the paper's prior (`nbsr-hmc`, `nbsr-laplace`), with the diffuse intercept prior (`*-wideint`) and the empirical-Bayes prior (`nbsr-laplace-eb`); LinDA, ALDEx2, ANCOM-BC, ANCOM-BC2; DESeq2 (with and without apeglm), edgeR QLF and LRT. Three files to stay under GitHub's file-size limit. |
| `microbiome_features.csv.gz` | microbiome, logistic-normal Dirichlet-multinomial (ALL25/ALL50/BAL50; n = 40, 80, 160, 320) | `nbsr-laplace-wideint` (the paper's model), the per-sample dispersion variants (`*-wsample_wideint`, `*-wdepth_wideint`), `nbsr-laplace-eb`, the HMC check (`nbsr-hmc-wideint`, ALL50 at n = 40 and 160, five replicates); LinDA, ALDEx2, ANCOM-BC2 |
| `microbiome_Ai_features.csv.gz` | the same generator with a concentration drawn per sample (ALL25, BAL50) | the three NBSR dispersion models |
| `conc_A<A>_features.csv.gz` | the concentration sweep: plain Dirichlet-multinomial, ALL25, A = 38, 100, 400, 1600, 4000 and A = 100 at ten times the depth | `nbsr-laplace-wideint` |
| `selex_features.csv.gz` | the selection experiment: `sample_7` is the full data, `sample_5` the 20 random 5 + 5 subsets | NBSR (Laplace, elicited prior `nbsr-laplace-eb`; fixed priors; HMC `nbsr-hmc-eb`) and the seven competitor configurations |
| `permutation_features.csv.gz` | label permutations of the Baxter (CRC) and Schubert (CDI) studies, 50 per n | `nbsr-laplace-eb` and the competitors |
| `real_labels_features.csv.gz` | the two studies with their real labels (BAXTER_CRC, BAXTER_ADENOMA, SCHUBERT_CDI) | `nbsr-laplace-eb` and the competitors |
| `ablation_oracle_features.csv.gz`, `ablation_priors_features.csv.gz`, `ablation_dispersion_features.csv.gz` | the miRNA ablations of the appendix: the oracle and offset dispersion models with and without the diffuse intercept prior, the reference-feature and no-pivot variants, the prior ablations (`sigma2`, `bpi1`, `nobj`, the horseshoe `nbsr_hs-hmc`, the group term `wtrt`), and the dispersion-structure variants (`nbsr_disp-hmc-*`) | NBSR variants only |

Each `*_manifest.csv` lists rows and units per method; each `*_units.csv` the
wall time (`wall_seconds`) and cores (`chains`) of every fit.

## Columns

Identity: `method`, `simul_case`, `sample_size` (per condition), `replicate_id`,
`feature`; microbiome rows also carry `stratum` (HIGH/MED/LOW by template mean
proportion) and the template's `mean_prop_template`, `prevalence_template`.

Truth (from the generator's `fc.csv`): `alpha_null`, `alpha_alt` are the
Dirichlet parameters of the two conditions, i.e. the absolute abundances;
`true_effect_abs = log(alpha_alt / alpha_null)` is the absolute-scale log fold
change, zero for every unperturbed feature, and `is_de` flags the perturbed
ones; `true_effect_clr` is that effect centred over the features the method
tested; `true_effect` is the target the method is scored against
(`target_scale`: `absolute` for every method but ALDEx2, `clr` for ALDEx2);
`log_fc` (`log2_fc`) is the change in *proportion*, which is non-zero for
unperturbed features and is not a scoring target. SELEX rows carry `nixon27`
(the 27 characterised variants), `ltr_survival`, `onu_survival`,
`minutes_to_midlog` instead of a simulated truth, and `true_effect_abs` is
`1` for the 27 as a flag only.

NBSR rows: `post_mean`, `post_sd`, `p_pos` (posterior probability that the
shift-adjusted effect is positive; `lfsr = min(p_pos, 1 - p_pos)`), the
posterior quantiles `q005` ... `q995`, `centring_fraction` (the reference set:
0.1 = most abundant tenth, floor 30, on the miRNA benchmark; 1 = whole panel
elsewhere), and the delta-method statistic where it was computed: `dm_est`,
`dm_se`, `p_dm`. The posterior-z Wald p is `2 * pnorm(-abs(post_mean /
post_sd))`, adjusted by Benjamini-Hochberg within a replicate.

Competitor rows: `estimate`, `se` (natural-log scale; log2 outputs converted;
edgeR's `se` is derived from its statistic and gives no interval), `p_value`,
`q_value` (the method's own adjusted p), `called` (its own rule), `passed_sens`
(ANCOM-BC2's sensitivity filter), `tested` (FALSE for features a method
declined: ANCOM-BC's prevalence cut, DESeq2's independent filtering),
`native_tested`; `aldex_effect`, `aldex_diff_btw` are ALDEx2's own effect
columns.

## This bundle

`paper_tables.R` and `make_paper_figures.R` rebuild the paper's tables and figures from these files
(`Rscript paper_tables.R <this directory> <output directory>`, likewise for the figures; they are the
scripts of `recomb2026/` in <https://github.com/junseonghwan/nbsr-experiments>, which also holds the
generators of the simulated inputs and the Stan models). `paper_inputs.tar.gz`, deposited beside this
directory, holds the simulated and prepared inputs of every fit (`Y.csv`, `X.csv`, `fc.csv` per
replicate).

Staged on 2026-10-21 for download and
a Zenodo deposit. `make_paper_figures.R` is the paper repository's script
(`scripts/make_paper_figures.R`); it reads `data/nbsr_outputs/` relative to the
paper directory, so either place these files there or change `data_dir` at
its top. `SHA256SUMS` covers every file here; check with `sha256sum -c
SHA256SUMS`. The fits themselves were produced by the pipeline at
<https://github.com/junseonghwan/nbsr-benchmarks> (`scripts/export_paper_data.R`
wrote these tables from the per-unit `summary.rds` files).
