# Additional simulation workflows

- [SuSiEx and MultiSuSiE baseline benchmark](susiex_multisusie_baseline/): Benchmarks SuSiEx and MultiSuSiE against INFUSE and other methods in the baseline simulation setting using credible sets, PIPs, calibration, and functional enrichment.
- [PIP calibration](pip_calibration/): Evaluates PIP calibration for INFUSE and competing methods under both baseline annotation settings.
- [Ancestry-configuration classification](configuration_classification/): The scripts extract known-causal SNP configuration probabilities, calculate classification summaries, and evaluate retained credible sets containing exactly one causal SNP.
- [Annotation heterogeneity](annotation_heterogeneity/): Sensitivity simulations assessing whether INFUSE benefits from functional annotations when annotation values or enrichment effects differ across ancestries.
- [LD mismatch](ld_mismatch/): Additional simulations evaluating LD outlier detection under local LD discrepancies and simulated population mismatch affecting one or both ancestries.
- [Variant availability](variant_availability/): Comparisons of INFUSE, INFUSE without functional annotations (INFUSE−), SuSiEx, and MultiSuSiE while retaining variants available in at least one ancestry, alongside single-ancestry SuSiE.

Shared baseline definitions and plotting helpers are in `common/`. Set input and output paths in the relevant scripts or environment variables.
