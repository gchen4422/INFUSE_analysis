.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x

NATIVE_DIR <- file.path(.infuse_work, "INFUSE_variant_availability/simulation_native_assv")
MODULE_DIR <- dirname(NATIVE_DIR)

SOURCE_REGION_FILE <- file.path(Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project)), "simulation_locus_info.txt")
SOURCE_BFILE <- c(
  EUR = Sys.getenv("INFUSE_EUR_BFILE", "/path/to/genotypes/EUR"),
  AFR = Sys.getenv("INFUSE_AFR_BFILE", "/path/to/genotypes/AFR")
)
PLINK2_BIN <- Sys.getenv("INFUSE_PLINK2_BIN", unname(Sys.which("plink2")))
ANNOT_SOURCE_DIR <- Sys.getenv("INFUSE_BASELINE_ANNOTATIONS", "/path/to/baselineLF2.2.UKB")

ANCESTRIES <- c("EUR", "AFR")
N_LOCI <- 100L

ANNOT_SIGNAL_COLS <- c("H3K27ac_Hnisz_common", "DHS_Trynka_common",
                       "H3K4me3_Trynka_common", "Conserved_LindbladToh_common",
                       "non_synonymous_common")
ANNOT_NULL_COLS <- c("Coding_UCSC_common", "Coding_UCSC.flanking.500_common",
                     "CTCF_Hoffman_common", "Enhancer_Andersson_common",
                     "DHS_peaks_Trynka_common")
ANNOT_FIT_COLS <- c(ANNOT_SIGNAL_COLS, ANNOT_NULL_COLS)
ENRICHMENT_BETA <- 2

REGION_FILE <- file.path(NATIVE_DIR, "metadata", "regions.tsv")
PANEL_ROOT <- file.path(NATIVE_DIR, "panels_source_preserved")
WORK_ROOT <- file.path(NATIVE_DIR, ".work_source_preserved")
LOG_ROOT <- file.path(NATIVE_DIR, "logs_source_preserved")

stop_if_missing_sources <- function() {
  required <- c(SOURCE_REGION_FILE, PLINK2_BIN,
                paste0(rep(SOURCE_BFILE, each = 3), rep(c(".bed", ".bim", ".fam"), 2)))
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing read-only source file(s):\n", paste(missing, collapse = "\n"))
  invisible(TRUE)
}
