.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# One explicit baseline definition shared by calibration and configuration.
figure2_baseline_spec <- function(setting) {
  base <- .infuse_project
  stopifnot(setting %in% c("d10", "d187"))
  high <- setting == "d187"
  root <- if (high) file.path(base, "More_annotes") else base
  list(
    dimension = if (high) "d = 187" else "d = 10",
    d = if (high) 187L else 10L,
    result_root = root,
    result_subdir = file.path("result", if (high) "w2_modified_llk" else "w2_startover_baseline"),
    data_root = root,
    data_subdir = file.path("summary_data", "w2_startover"),
    figure_summary = file.path(root, "res_summary", if (high)
      "shared_50_more_annotes_mlk_round2_w2_startover_sparsepca_bic_modified_llk.RData"
      else "shared_50_sample_mlk_w2_startover.RData"),
    # Historical filename, not an assertion that d=10 used SPCA. Its per-SNP
    # data were checked against every raw baseline cell before use.
    legacy_cache = file.path(root, "res_summary", "shared_50_sample_mlk_w2_startover_spca.RData"),
    verified_cache = file.path(.infuse_work, "baseline_figure2_alignment_20260907/inputs", paste0(setting, "_figure2_verified_pip_data.rds"))
  )
}
