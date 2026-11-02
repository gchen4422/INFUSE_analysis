# INFUSE Analysis

Analysis code accompanying the **INFUSE** manuscript: INFUSE enables robust and powerful multi-ancestry fine-mapping by modeling ancestry-specific variants and integrating high-dimensional functional annotations.

This repository contains scripts for the simulation studies, LD-discrepancy analyses, and All of Us real-data applications. The core INFUSE method is implemented in the separate [`INFUSE`](https://github.com/gchen4422/INFUSE) software package.

## Repository Structure

```text
INFUSE_analysis/
├── Simulation/
│   ├── informative/            # Selected annotation setting: 10 annotations
│   ├── relevance_uncertain/    # High-dimensional setting: 187 baselineLF2.2.UKB annotations
│   ├── external_ld/            # External LD reference analyses using 1000 Genomes Phase 3
│   └── revision/               # Additional simulation benchmarks and sensitivity analyses
└── Real_data/                  # All of Us BMI, DBP, and SBP analyses
    ├── 5_run_finemap/
    ├── evo2/
    └── revision/               # Variant-availability and Evo 2 analyses
```

## Simulation

Simulations use UK Biobank imputed genotypes to generate multi-ancestry GWAS summary statistics across 100 genomic regions. Benchmark methods include INFUSE, MESuSiE, SuSiEx, MultiSuSiE, ancestry-specific SuSiE with post-processing, PAINTOR, and PAINTOR with annotations. Analyses include selected-informative and high-dimensional relevance-uncertain annotation settings, varying annotation enrichment strengths and cross-ancestry causal sharing levels, as well as robustness scenarios for sample imbalance, external 1000 Genomes Phase 3 LD reference panels, and simulated LD-mismatch/outlier variants.

### `Simulation/informative/`

Selected informative annotation setting with 10 annotations, including enriched and non-informative annotations.

| Script                                  | Description                                                   |
| --------------------------------------- | ------------------------------------------------------------- |
| `1_sample_locus.R`                      | Sample candidate genomic loci                                 |
| `2_locus_selection.R`                   | Define fine-mapping regions                                   |
| `3_shared_50_simulation_mesusie_inf*.R` | Simulate multi-ancestry GWAS summary statistics               |
| `4_shared_50_*.R/.sh`                   | Run fine-mapping methods                                      |
| `5_result_summarizing*.R`               | Summarize power, FDR, credible sets, calibration, and runtime |
| `6_enrichment_prepration.R`             | Prepare enrichment inputs                                     |
| `7_enrichment_results_summary.R`        | Summarize enrichment results                                  |
| `8_*enrichmen_results*.R`               | Generate enrichment plots                                     |
| `utility*.R`                            | Shared helper functions                                       |

### `Simulation/relevance_uncertain/`

High-dimensional annotation setting using 187 baselineLF2.2.UKB annotations. This folder also includes sparse-PCA annotation processing and LD-discrepancy analyses.


| Script / Folder                                  | Description                                                                        |
| ------------------------------------------------ | ---------------------------------------------------------------------------------- |
| `3_shared_50_simulation_mesusie_inf_20_annotes/` | Simulate high-dimensional summary statistics and LD-mismatch scenarios             |
| `4_shared_50_mesusie_all_annotes_runfinemap/`    | Run fine-mapping methods                                                           |
| `5_result_summarizing/`                          | Summarize power, FDR, credible sets, calibration, runtime, and LD-mismatch results |
| `6_enrichment_prepration.R`                      | Prepare enrichment inputs                                                          |
| `7_enrichment_results_summary.R`                 | Summarize enrichment results                                                       |
| `8_estimate_toppcs_nonzeros.R`                   | Explore sparse PCs and non-zero annotation loadings                                |
| `9_*enrichmen_results*.R`                        | Generate enrichment plots                                                          |
| `10_test_ld_outlier/`                            | Evaluate the multi-ancestry LD-discrepancy test                                    |
| `Annotation_selection.R`                         | Process annotation inputs                                                          |
| `Categorize_baseline_annot.R`                    | Categorize baselineLF2.2.UKB annotations                                           |


### `Simulation/external_ld/`

External LD reference analyses using 1000 Genomes Phase 3 LD panels.

| Script                                 | Description                 |
| -------------------------------------- | --------------------------- |
| `3_shared_50_simulation_mesusie_inf.R` | Simulate summary statistics |
| `4_shared_50_*.R/.sh`                  | Run fine-mapping methods    |
| `5_result_summarizing.R`               | Summarize results           |
| `utility.R`                            | Shared helper functions     |

### `Simulation/revision/`

Additional benchmarks and sensitivity analyses covering SuSiEx and MultiSuSiE, PIP calibration, ancestry-configuration classification, annotation heterogeneity, LD mismatch, and ancestry-specific variant availability. See the [simulation revision README](Simulation/revision/README.md) for workflow details and running instructions.

## Real Data: All of Us BMI, DBP, and SBP

Real-data analyses use All of Us v8 GWAS summary statistics for BMI, DBP, and SBP in European- and African-ancestry cohorts. Methods include INFUSE, MESuSiE, SuSiEx, PAINTOR, and PAINTOR with annotations. Evo 2-derived variant scores are evaluated as additional annotations.

| Script                           | Description                                                 |
| -------------------------------- | ----------------------------------------------------------- |
| `1_prepare_for_finemap.R`        | Format and QC GWAS summary statistics                       |
| `2_define_region_1mb.R`          | Define candidate regions                                    |
| `3_generate_loci.R`              | Generate fine-mapping loci                                  |
| `4_generate_loci_sumstats.R/.sh` | Extract per-locus summary statistics and LD inputs          |
| `5_run_finemap/*`                | Run fine-mapping methods by trait                           |
| `6_annotation.R`                 | Annotate prioritized variants                               |
| `7_result_summarizing*.R`        | Summarize fine-mapping and Evo 2 results                    |
| `8_plot_stat*.R`                 | Generate summary plots                                      |
| `9_plot_enrichment*.R`           | Plot enrichment for variants prioritized in either ancestry |
| `10_plot_enrichment*.R`          | Plot shared and ancestry-specific enrichment                |
| `11_plot_locus_zoom.R`           | Generate locus-level visualizations                         |
| `evo2/score_evo2_*.sh`           | Compute Evo 2-7B variant scores                             |
| `utility.R`                      | Shared helper functions                                     |
| [`revision/`](Real_data/revision/README.md) | Variant-availability and Evo 2 analyses; see its README for details |

## Dependencies

**R packages:** `INFUSE`, `MESuSiE`, `susieR`, `data.table`, `dplyr`, `ggplot2`, `patchwork`, `mvtnorm`, `forcats`, `GenomicScores`, `Matrix`, `reticulate`, `snpStats`, `tidyr`, `ggpattern`, `ggrepel`, `gridExtra`, `egg`, `mixsqp`

**Python packages:** `MultiSuSiE`, `numpy`, `pandas`, `matplotlib`, `pyarrow`

**External tools:** SuSiEx, PAINTOR, PLINK/PLINK2, REGENIE, Evo 2-7B

**Required data:** UK Biobank genotypes, All of Us v8 GWAS data, 1000 Genomes Phase 3 reference panels, baselineLF2.2.UKB annotations, hg38 reference genome, and Evo 2-7B resources.

Controlled-access datasets and large reference files are not included in this repository.

## Notes

Scripts assume a SLURM-based HPC environment. File paths are cluster-specific and may need to be updated before reuse. This repository documents the manuscript analysis workflow and is not intended as a fully containerized one-command reproduction pipeline.
