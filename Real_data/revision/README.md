# Additional real-data workflows

- [Variant availability](variant_availability/): All of Us fine-mapping retaining variants available in at least one ancestry.
- [Evo 2](evo2/): Evaluates predicted functional disruption among variants prioritized by INFUSE and MESuSiE using raw signed Evo 2 scores.

Code for the REVEL/SpliceAI external-annotation analyses is available upon request for use within the All of Us Researcher Workbench.

Set input and output paths in the relevant scripts or environment variables. Controlled-access inputs must be supplied separately.

## Variant availability

This workflow extends the All of Us BMI, DBP, and SBP analyses to retain variants available in at least one ancestry. It compares INFUSE, SuSiEx, and INFUSE without functional annotations (INFUSE−) using credible-set resolution and functional enrichment. Ancestry-specific SuSiE is also included as a comparator.

### Workflow order

The main scripts in `variant_availability/` are numbered continuously in execution order:

1. `01_extract_baseline.py`
2. `02_match_parquet.py`
3. `03_reuse_evaluation_annotations.R`
4. `04_freeze_variant_manifest.R`
5. `05_prepare_locus.R`
6. `06_run_method.R`
7. `07_run_susiex.R`
8. `08_summarize_results.R`
9. `09_matched_size_enrichment.R`
10. `10_pooled_configuration_enrichment.R`
11. `11_plot_enrichment.R`
12. `12_finalize_matched_summary.R`
13. `13_cs_summary_panels.R`
14. `14_combine_real_data_panels.R`
15. `15_final_bmi_bp_figure.R`

`variant_availability/run_summary.sh` runs the summary/plotting stages in their required order;
`14_combine_real_data_panels.R` also invokes the final figure script.
`run_job.sh`, `run_array.sh`, `run_matched_size.sh`, `plink_limited.sh` and
`annotate_enrichment_panels.R` are unnumbered execution/annotation helpers.
The launchers call the source files in `variant_availability/` directly. Supply the existing post-QC manifests and observed-input objects in the chosen All of Us output/input tree before starting; the controlled-access input-audit pipeline is not distributed.

## Evo 2 external validation

`evo2/` compares raw signed Evo 2 scores among prioritized and background variants, with locus-cluster bootstrap confidence intervals and summary plots. More negative scores indicate greater predicted functional disruption. Embeddings and raw score datasets are not included.

### Scripts

1. `01_add_raw_evo2_scores.R`
2. `02_plot_enrichment_either_bp_evo2.R`
3. `03_plot_enrichment_either_bmi_evo2.R`
4. `04_plot_bmi_asinh_sensitivity.R`
5. `05_export_allofus_annotation_benchmark.R`
6. `06_create_primary_supp_figure.R`
7. `plot_evo2_with_ci.py`

`00_evo2_enrichment_common.R` provides shared helpers; `submit_*.sh` scripts submit the corresponding analyses.
