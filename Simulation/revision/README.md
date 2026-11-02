# Additional simulation workflows

- [SuSiEx and MultiSuSiE baseline benchmark](susiex_multisusie_baseline/): Runs the additional SuSiEx and MultiSuSiE fits in the baseline simulation setting, followed by credible-set, PIP, calibration, and enrichment summaries.
- [PIP calibration](pip_calibration/): Summarizes baseline PIP calibration for both annotation dimensions.
- [Ancestry-configuration classification](configuration_classification/): The scripts extract known-causal SNP configuration probabilities, calculate classification summaries, and evaluate retained credible sets containing exactly one causal SNP.
- [Annotation heterogeneity](annotation_heterogeneity/): Evaluates performance when annotation values and their effects on causal probability differ across ancestries.
- [LD mismatch](ld_mismatch/): Evaluates LD outlier detection under local LD discrepancies and simulated population mismatch affecting one or both ancestries.
- [Variant availability](variant_availability/): Comparisons of INFUSE, INFUSE without functional annotations (INFUSE−), SuSiEx, and MultiSuSiE while retaining variants available in at least one ancestry, alongside single-ancestry SuSiE.

Shared baseline definitions and plotting helpers are in `common/`. Set input and output paths in the relevant scripts or environment variables.
