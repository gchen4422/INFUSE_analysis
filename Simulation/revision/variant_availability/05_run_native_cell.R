#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) < 4L) {
  stop("Usage: Rscript 05_run_native_cell.R <locus> <base_index> <h2_index> <n_assv_each> [L] [run_susiex]")
}
locus <- as.integer(args[1L]); base_index <- as.integer(args[2L])
h2_index <- as.integer(args[3L]); n_assv_each <- as.integer(args[4L])
L_effects <- if (length(args) >= 5L) as.integer(args[5L]) else 10L
run_sx <- if (length(args) >= 6L) as.logical(as.integer(args[6L])) else TRUE

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                    value = TRUE)[1L]), mustWork = TRUE)
NATIVE_SRC_DIR <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_sim_common.R"))
suppressMessages({
  library(INFUSE)
  library(INFUSE)
  library(MESuSiE)
  library(susieR)
})

stopifnot(locus %in% 1:100, base_index %in% 1:3, h2_index %in% 1:2,
          n_assv_each %in% NATIVE_ASSV_LEVELS)
out_file <- native_cell_path(locus, base_index, h2_index, n_assv_each)
if (file.exists(out_file)) {
  cat("complete output exists; leaving untouched:", out_file, "\n")
  quit(status = 0L)
}
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(NATIVE_SIM_WORK, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("native cell locus=%d base=%d h2=%g assv_each=%d L=%d SuSiEx=%s\n",
            locus, NATIVE_BASE_COUNTS[base_index], NATIVE_H2[h2_index],
            n_assv_each, L_effects, run_sx))
lo <- load_native_locus(locus, load_genotypes = TRUE)
design <- draw_native_truth(lo, base_index)
sim <- simulate_native_z(lo, design, base_index, h2_index, n_assv_each)

# Structural assertions: native AS-SVs exist only in their causal ancestry.
tr <- sim$truth
stopifnot(sum(tr$causal_origin == "native_ASSV_EUR") == n_assv_each,
          sum(tr$causal_origin == "native_ASSV_AFR") == n_assv_each,
          all(tr$available_EUR[tr$causal_origin == "native_ASSV_EUR"]),
          all(!tr$available_AFR[tr$causal_origin == "native_ASSV_EUR"]),
          all(!tr$available_EUR[tr$causal_origin == "native_ASSV_AFR"]),
          all(tr$available_AFR[tr$causal_origin == "native_ASSV_AFR"]),
          sum(tr$Signal != 0L) == NATIVE_BASE_COUNTS[base_index] + 2L * n_assv_each)

sx_work <- file.path(NATIVE_SIM_WORK,
                     sprintf("locus_%03d_base%d_h2%d_assv%d_%d", locus,
                             base_index, h2_index, n_assv_each, Sys.getpid()))
mm <- run_methods_v2(sim$zfile, lo$avail, lo$R_union, lo$annot, L = L_effects,
                     map = lo$map, locus = locus, susiex_work = sx_work,
                     run_susiex_flag = run_sx)

expected <- c("infuse_union", "mesusie_union", "infuse_intersect",
              "mesusie_intersect", "susie")
if (run_sx) expected <- c(expected, "susiex_intersect", "susiex_union")
missing_methods <- setdiff(expected, names(mm$res))
if (length(missing_methods)) {
  stop("Method failure(s): ", paste(missing_methods, collapse = ", "),
       "; scratch preserved at ", sx_work)
}
unlink(sx_work, recursive = TRUE)

# Intersection methods must assign exact zero to every native AS-SV causal.
native_idx <- which(grepl("^native_ASSV", tr$causal_origin))
for (m in c("infuse_intersect", "mesusie_intersect", "susiex_intersect")) {
  if (m %in% names(mm$res) && length(native_idx) && any(mm$res[[m]]$pip[native_idx] != 0))
    stop(m, " assigned nonzero PIP to a native AS-SV")
}

cell <- list(
  metadata = list(locus = locus, base_count_index = base_index,
                  n_base = NATIVE_BASE_COUNTS[base_index], h2_index = h2_index,
                  h2_per_causal = NATIVE_H2[h2_index], n_assv_each = n_assv_each,
                  L = L_effects, run_susiex = run_sx),
  truth = tr, availability = lo$avail, zfile = sim$zfile,
  causal_beta = sim$beta, h2_realized = sim$h2_realized,
  methods = mm$res, timing_minutes = mm$timing)
tmp <- paste0(out_file, ".partial_", Sys.getpid())
saveRDS(cell, tmp, compress = FALSE)
if (!file.rename(tmp, out_file)) stop("Atomic output rename failed: ", out_file)
cat("wrote", out_file, "\n")
