#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("d10", "d187")) {
  stop("Usage: 01_generate_calibration.R <d10|d187>", call. = FALSE)
}
setting <- args[[1]]

base <- .infuse_project
revision <- .infuse_work
analysis_dir <- file.path(revision, "pip_calibration_w2_d10_d187")
figure_root <- file.path(revision, "figures", "pip_calibration_w2_d10_d187")
source(file.path(.infuse_repo, "Simulation/revision/common/baseline_sources.R"))
baseline_spec <- figure2_baseline_spec(setting)

settings <- list(
  d10 = list(
    label = "d10_direct_mlk",
    d = 10L,
    fit = "direct MLK (no sparse PCA)",
    input = figure2_baseline_spec("d10")$verified_cache
  ),
  d187 = list(
    label = "d187_sparse_pca",
    d = 187L,
    fit = "Figure 2 high-dimensional MLK baseline (w2_modified_llk)",
    input = figure2_baseline_spec("d187")$verified_cache
  )
)
cfg <- settings[[setting]]
out_dir <- file.path(figure_root, cfg$label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(data.table)
})
source(file.path(.infuse_repo, "Simulation/informative/utility.R"))

verified <- readRDS(cfg$input)
if (!identical(verified$baseline$result_subdir, baseline_spec$result_subdir) ||
    !identical(verified$baseline$result_root, baseline_spec$result_root) ||
    verified$verified_cells != 600L || is.null(verified$data_all)) {
  stop("Input does not match the Figure 2 baseline: ", cfg$input, call. = FALSE)
}
data_all <- verified$data_all
rm(verified)
gc()

required <- c(
  "Signal", "h2", "causal_num", "MESuSiE_Either", "SuSiE_Either",
  "Paintor_PIP", "MAFPC_Either", "Paintor_annot_pip", "MESuSiE_Shared",
  "SuSiE_Shared", "MAFPC_Shared", "MESuSiE_WB", "MESuSiE_BB",
  "SuSiE_WB", "SuSiE_BB", "MAFPC_WB", "MAFPC_BB"
)
missing <- setdiff(required, names(data_all))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
stopifnot(nrow(data_all) == 1314114L,
          identical(as.integer(sort(unique(data_all$Signal))), 0:3),
          all(data_all$h2 %in% 1:2),
          all(data_all$causal_num %in% 1:3))

relabel_methods <- function(x, mapping, order) {
  old <- as.character(x)
  new <- unname(mapping[old])
  if (anyNA(new)) stop("Unmapped method label(s): ", paste(unique(old[is.na(new)]), collapse = ", "))
  factor(new, levels = order)
}

# Either configuration: causal in EUR, AFR, or both.
either_cols <- c("MESuSiE_Either", "SuSiE_Either", "Paintor_PIP",
                 "MAFPC_Either", "Paintor_annot_pip")
cal_either <- create_obs_frq_byh2(
  data_all %>% select(Signal, h2, causal_num, all_of(either_cols)),
  c(1, 2, 3), either_cols
)
either_map <- c(
  MESuSiE_Either = "MESuSiE", SuSiE_Either = "SuSiE",
  Paintor_PIP = "Paintor", MAFPC_Either = "INFUSE",
  Paintor_annot_pip = "Paintor+annot"
)
method_order <- c("MESuSiE", "SuSiE", "Paintor", "INFUSE", "Paintor+annot")
cal_either$Method <- relabel_methods(cal_either$Method, either_map, method_order)

# Shared configuration: only Signal == 3 is a positive. The legacy source used
# c(1,2,3), which calibrates a shared-configuration PIP against the wrong truth.
shared_cols <- c("MESuSiE_Shared", "SuSiE_Shared", "Paintor_PIP",
                 "MAFPC_Shared", "Paintor_annot_pip")
cal_shared <- create_obs_frq_byh2(
  data_all %>% select(Signal, h2, causal_num, all_of(shared_cols)),
  3, shared_cols
)
shared_map <- c(
  MESuSiE_Shared = "MESuSiE", SuSiE_Shared = "SuSiE",
  Paintor_PIP = "Paintor", MAFPC_Shared = "INFUSE",
  Paintor_annot_pip = "Paintor+annot"
)
cal_shared$Method <- relabel_methods(cal_shared$Method, shared_map, method_order)

# Ancestry-specific configurations: Signal == 1 is WB/EUR and Signal == 2 is
# BB/AFR. The legacy snippet accidentally used Signal == 1 for both panels.
data_all <- data_all %>% mutate(
  Paintor_WB = Paintor_PIP,
  Paintor_BB = Paintor_PIP,
  Paintor_annot_WB = Paintor_annot_pip,
  Paintor_annot_BB = Paintor_annot_pip
)
wb_cols <- c("MESuSiE_WB", "Paintor_WB", "SuSiE_WB", "MAFPC_WB", "Paintor_annot_WB")
bb_cols <- c("MESuSiE_BB", "Paintor_BB", "SuSiE_BB", "MAFPC_BB", "Paintor_annot_BB")
# Match either/shared: pool SNP records across causal counts within each h2.
cal_wb <- create_obs_frq_byh2(
  data_all %>% select(Signal, h2, causal_num, all_of(wb_cols)),
  1, wb_cols
)
cal_bb <- create_obs_frq_byh2(
  data_all %>% select(Signal, h2, causal_num, all_of(bb_cols)),
  2, bb_cols
)
cal_ancestry <- bind_rows(cal_wb, cal_bb)

ancestry_map <- c(
  MESuSiE_WB = "MESuSiE~White~British", MESuSiE_BB = "MESuSiE~Black~British",
  Paintor_WB = "Paintor~White~British", Paintor_BB = "Paintor~Black~British",
  SuSiE_WB = "SuSiE~White~British", SuSiE_BB = "SuSiE~Black~British",
  MAFPC_WB = "INFUSE~White~British", MAFPC_BB = "INFUSE~Black~British",
  Paintor_annot_WB = "Paintor+annot~White~British",
  Paintor_annot_BB = "Paintor+annot~Black~British"
)
ancestry_order <- unname(ancestry_map[c(
  "MESuSiE_WB", "MESuSiE_BB", "Paintor_WB", "Paintor_BB", "SuSiE_WB",
  "SuSiE_BB", "MAFPC_WB", "MAFPC_BB", "Paintor_annot_WB", "Paintor_annot_BB"
)])
cal_ancestry$Method <- relabel_methods(cal_ancestry$Method, ancestry_map, ancestry_order)

p_either <- PIP_calibration_shared_byh2_fun(cal_either)
p_shared <- PIP_calibration_shared_byh2_fun(cal_shared)
p_ancestry <- PIP_calibration_ancestry_byh2_fun(cal_ancestry) +
  theme(strip.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "pip_either_calibration.pdf"), p_either,
       height = 200, width = 300, units = "mm", dpi = 500)
ggsave(file.path(out_dir, "pip_shared_calibration.pdf"), p_shared,
       height = 200, width = 300, units = "mm", dpi = 500)
ggsave(file.path(out_dir, "pip_ancestry_calibration.pdf"), p_ancestry,
       height = 200, width = 400, units = "mm", dpi = 500)

fwrite(as.data.frame(cal_either), file.path(out_dir, "pip_either_calibration.tsv"), sep = "\t")
fwrite(as.data.frame(cal_shared), file.path(out_dir, "pip_shared_calibration.tsv"), sep = "\t")
fwrite(as.data.frame(cal_ancestry), file.path(out_dir, "pip_ancestry_calibration.tsv"), sep = "\t")


