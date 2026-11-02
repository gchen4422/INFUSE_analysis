#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

revision <- .infuse_work
figure_root <- file.path(revision, "figures", "pip_calibration_w2_d10_d187")
d10_dir <- file.path(figure_root, "d10_direct_mlk")
d187_dir <- file.path(figure_root, "d187_sparse_pca")
out_dir <- file.path(figure_root, "primary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c(
  "SuSiE", "Paintor", "Paintor+", "Paintor+H",
  "MESuSiE", "INFUSE", "INFUSE-H"
)

read_endpoint <- function(setting_dir, endpoint) {
  path <- file.path(setting_dir, sprintf("pip_%s_calibration.tsv", endpoint))
  if (!file.exists(path)) stop("Missing calibration table: ", path, call. = FALSE)
  fread(path)
}

combine_nonancestry <- function(endpoint) {
  d10 <- read_endpoint(d10_dir, endpoint)
  d187 <- read_endpoint(d187_dir, endpoint)

  # Baseline methods are represented once, using the d=10 calibration table.
  baseline <- d10[Method %chin% c("SuSiE", "Paintor", "MESuSiE")]
  baseline[, SourceSetting := "d10 baseline"]

  informative <- d10[Method %chin% c("Paintor+annot", "INFUSE")]
  informative[, Method := fifelse(Method == "Paintor+annot", "Paintor+", "INFUSE")]
  informative[, SourceSetting := "d10 informative annotations"]

  high_dimensional <- d187[Method %chin% c("Paintor+annot", "INFUSE")]
  high_dimensional[, Method := fifelse(Method == "Paintor+annot", "Paintor+H", "INFUSE-H")]
  high_dimensional[, SourceSetting := "d187 high-dimensional annotations"]

  out <- rbindlist(list(baseline, informative, high_dimensional), use.names = TRUE)
  out[, Method := factor(Method, levels = method_levels)]
  if (anyNA(out$Method) || !identical(levels(out$Method), method_levels)) {
    stop("Failed to construct the requested method ordering.", call. = FALSE)
  }
  out[, h2 := factor(
    as.character(h2),
    levels = c("~h^2 == 10^-4", "~h^2 == 2%*%10^-4"),
    labels = c("h^2 == 10^-4", "h^2 == 2 %*% 10^-4")
  )]
  if (anyNA(out$h2)) stop("Unexpected h2 label.", call. = FALSE)
  setorder(out, Method, h2, gr)
  out[]
}

parse_ancestry <- function(x) {
  x <- as.character(x)
  ancestry <- fifelse(
    grepl("~White~British$", x), "White British",
    fifelse(grepl("~Black~British$", x), "Black British", NA_character_)
  )
  base_method <- sub("~(White|Black)~British$", "", x)
  data.table(BaseMethod = base_method, Ancestry = ancestry)
}

combine_ancestry <- function() {
  d10 <- read_endpoint(d10_dir, "ancestry")
  d187 <- read_endpoint(d187_dir, "ancestry")
  d10 <- cbind(d10, parse_ancestry(d10$Method))
  d187 <- cbind(d187, parse_ancestry(d187$Method))
  if (anyNA(d10$Ancestry) || anyNA(d187$Ancestry)) {
    stop("Could not parse an ancestry method label.", call. = FALSE)
  }

  baseline <- d10[BaseMethod %chin% c("SuSiE", "Paintor", "MESuSiE")]
  baseline[, `:=`(Method = BaseMethod, SourceSetting = "d10 baseline")]

  informative <- d10[BaseMethod %chin% c("Paintor+annot", "INFUSE")]
  informative[, Method := fifelse(BaseMethod == "Paintor+annot", "Paintor+", "INFUSE")]
  informative[, SourceSetting := "d10 informative annotations"]

  high_dimensional <- d187[BaseMethod %chin% c("Paintor+annot", "INFUSE")]
  high_dimensional[, Method := fifelse(BaseMethod == "Paintor+annot", "Paintor+H", "INFUSE-H")]
  high_dimensional[, SourceSetting := "d187 high-dimensional annotations"]

  out <- rbindlist(list(baseline, informative, high_dimensional), use.names = TRUE)
  out[, Method := factor(Method, levels = method_levels)]
  if (anyNA(out$Method)) stop("Unexpected ancestry method label.", call. = FALSE)

  out[, h2 := factor(
    as.character(h2),
    levels = c("~h^2 == 10^-4", "~h^2 == 2%*%10^-4"),
    labels = c("h^2 == 10^-4", "h^2 == 2 %*% 10^-4")
  )]
  out[, Ancestry := factor(Ancestry, levels = c("White British", "Black British"))]
  if (anyNA(out$h2) || anyNA(out$Ancestry)) {
    stop("Failed to construct h2-by-ancestry panel rows.", call. = FALSE)
  }
  stopifnot(!"causal_num" %in% names(out),
            nrow(unique(out[, .(h2, Ancestry)])) == 4L)
  setorder(out, Method, h2, Ancestry, gr)
  out[]
}

calibration_theme <- function() {
  theme_bw(base_size = 11) +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 9),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = grid::unit(3, "mm"),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      plot.caption = element_text(size = 8, hjust = 0)
    )
}

annotation_caption <- paste0(
  "'+' denotes the informative-annotation setting (d = 10); ",
  "'H' denotes the high-dimensional relevance-uncertain setting (d = 187)."
)

base_calibration_layers <- list(
  geom_pointrange(
    aes(ymin = obs_min, ymax = obs_max),
    colour = "grey55", linewidth = 0.3
  ),
  geom_point(colour = "black", size = 1.5),
  geom_abline(intercept = 0, slope = 1, colour = "red", linewidth = 0.4),
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ),
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ),
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)),
  labs(x = "Mean PIP", y = "Observed Frequency", caption = annotation_caption),
  calibration_theme()
)

make_nonancestry_plot <- function(dat) {
  ggplot(dat, aes(x = mean_pip, y = obs_frq)) +
    base_calibration_layers +
    facet_grid(
      rows = vars(h2), cols = vars(Method),
      labeller = labeller(h2 = label_parsed, Method = label_value)
    )
}

make_ancestry_plot <- function(dat) {
  ggplot(dat, aes(x = mean_pip, y = obs_frq)) +
    base_calibration_layers +
    facet_grid(
      rows = vars(h2, Ancestry), cols = vars(Method),
      labeller = labeller(h2 = label_parsed, Ancestry = label_value, Method = label_value)
    ) +
    labs(caption = paste0(annotation_caption,
      "\nWithin each h² and ancestry, SNP records are pooled across 2, 3, and 5 causal SNP settings."))
}

either <- combine_nonancestry("either")
shared <- combine_nonancestry("shared")
ancestry <- combine_ancestry()

expected_nonancestry <- method_levels
stopifnot(
  identical(levels(either$Method), expected_nonancestry),
  identical(levels(shared$Method), expected_nonancestry),
  identical(levels(ancestry$Method), expected_nonancestry)
)

p_either <- make_nonancestry_plot(either)
p_shared <- make_nonancestry_plot(shared)
p_ancestry <- make_ancestry_plot(ancestry)

ggsave(file.path(out_dir, "pip_either_calibration.pdf"), p_either,
       width = 360, height = 205, units = "mm", dpi = 500)
ggsave(file.path(out_dir, "pip_either_calibration.png"), p_either,
       width = 360, height = 205, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "pip_shared_calibration.pdf"), p_shared,
       width = 360, height = 205, units = "mm", dpi = 500)
ggsave(file.path(out_dir, "pip_shared_calibration.png"), p_shared,
       width = 360, height = 205, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "pip_ancestry_calibration.pdf"), p_ancestry,
       width = 360, height = 260, units = "mm", dpi = 500)
ggsave(file.path(out_dir, "pip_ancestry_calibration.png"), p_ancestry,
       width = 360, height = 260, units = "mm", dpi = 300)

fwrite(as.data.frame(either), file.path(out_dir, "pip_either_calibration.tsv"), sep = "\t")
fwrite(as.data.frame(shared), file.path(out_dir, "pip_shared_calibration.tsv"), sep = "\t")
fwrite(as.data.frame(ancestry), file.path(out_dir, "pip_ancestry_calibration.tsv"), sep = "\t")


