# Re-analysis and simulation for the authors' reply (November 2026)

Code for the re-analyses and simulations reported in the authors' reply to correspondence on

> Baumeister S-E, McGrath C, Kocher T, Nascimento GG. Childhood famine exposure and late-life
> edentulism: evidence from natural experiments. *Journal of Clinical Periodontology* 2026;53:1673–80.
> doi:10.1111/jcpe.70190

The scripts extend the published pipeline (`../famine_edentulism_analysis.R`). They address three
questions: (1) how graded exposure around the regression-discontinuity cut-off and susceptibility
beyond age 5 affect each design; (2) which effect sizes the designs, alone and combined, can exclude;
and (3) whether the results depend on the CLHLS outcome definition, on pooling famines of different
severity, or on provincial famine intensity within China.

**No individual-level data are included.** The scripts read the pooled analysis file written by the
published pipeline and, for `00_prep_data.R`, the raw survey files (obtain them from the providers;
see the main `README.md`). Local extracts (`data_local/`) and all outputs (`output/`) are excluded by
`.gitignore`.

## Running

Run every script from this folder, in numerical order (R ≥ 4.4):

```
Rscript 00_prep_data.R      # needs FAMINE_DATA_DIR (raw survey files)
Rscript 02_reproduce.R
...
```

Paths: `01_common.R` finds the repository root as the parent folder (or set `FAMINE_REPO_DIR`); the
pooled file `famine_edent_pooled.rds` must have been written there by the published pipeline. Packages:
dplyr, splines, sandwich, lmtest, haven, fixest, rdrobust (4.0.0), RDHonest (1.0.1), metafor, lme4,
parallel.

## Scripts

| Script | Content |
|---|---|
| `00_prep_data.R` | Local extracts: CLHLS tooth counts and birth province; CHARLS birth province |
| `01_common.R` | Estimators as precomputed closures (RDiT local-linear with covariates, reduced form, ITS variants, dose-response by famine years at ages 1–5); paths |
| `02_reproduce.R` | Reproduces the published RDiT (+0.134 pp) and reduced form (+0.824 pp); dose-response analysis |
| `03_its_specifications.R` | Interrupted time series: published specification and alternatives (cell weights, trend, window); lack-of-fit test |
| `04_plasmode.R`, `04b_gls.R` | Plasmode simulation (no-effect data-generating process fitted to the real data; injected uniform, graded and recovery-spill-over effects); GLS combination of the designs |
| `04c_plasmode_susceptibility.R` | Plasmode with susceptibility extending beyond age 5 |
| `05_outcome_threshold.R` | CLHLS edentulism coded as 0, ≤ 2 or ≤ 4 remaining teeth; harmonised-outcome subset (nine surveys) |
| `06_severity_did.R` | Region-of-birth famine-intensity difference-in-differences (CHARLS, CLHLS) with cluster-robust, individual-level and placebo-province permutation inference |
| `07_heterogeneity_severity.R` | Study-specific estimates, random-effects meta-analysis and meta-regression (minimum intake, duration, famine-to-survey interval, outcome type), duration strata, leave-one-out |
| `09_clhls_famine_years.R` | CLHLS famine dated 1942–43 instead of 1942–45 |
| `10_rd_inference.R` | rdrobust conventional, bias-corrected and robust inference; RDHonest honest intervals over a grid of curvature bounds |
| `11_severity_did_continuous.R` | Graded severity difference-in-differences in CHARLS by 1960 provincial excess mortality; permutation inference; event study |
| `12_extra_checks.R` | Benchmarks at the CHARLS baseline; information share by study; event-study trend tests |
| `13_verification_followups.R` | Equivalence tests; permutation-based intervals; discontinuity at the famine-end boundary; relative scale; survivor-selection bias grid; CLHLS record checks |
| `14_schoolage_severity_did.R` | Graded severity difference-in-differences for school-age exposure |
| `15_final_checks.R` | Honest RDiT on an adjusted outcome; cohort-level inference for the severity difference-in-differences; agreement test and weights of the GLS combination; duration strata with CLHLS as a short famine |
| `16_gls_sensitivity.R` | GLS combination leaving out each design and with inflated variances |
| `17_plasmode_rd_overnight.R`, `17b_plasmode_rd_summary.R` | Second simulation run (10 cells × 288 replications; about 14 hours on 12 cores): narrow-bandwidth conventional RDiT, rdrobust with the bandwidth re-selected in every replication, two honest fits; cohort-shock, relative-effect and non-Chinese-famine scenarios; resumable blocks |
| `18_sensitivity_reanalysis.R`, `18h_its_placebo_cutoffs.R` | Placebo cut-offs and placebo famine dates (calibration of the standard errors on the real data), ITS lack of fit, dose windows, leave-one-province-out, honest M grid, older-age boundaries |

## Input

`inputs/linyang2000_table3_province_death_rates.csv`: provincial crude death rates 1954–1966 as printed
in Lin JY, Yang DT (2000) Food availability, entitlements and the Chinese famine of 1959–61. *Economic
Journal* 110:136–158 (Table 3), doi:10.1111/1468-0297.00494, with derived excess death rates. The
province labels are kept as printed (`Shansi` = Shanxi, printed `Shanxi` = Shaanxi); the scripts map
them explicitly and check the mapping.

## Seeds

20261109 (04, 06, 11), 20261110 (04c), 20261111 (13), 20261112 (14), 20261113 (15),
20261114 + 7919 × cell + replication (17).
