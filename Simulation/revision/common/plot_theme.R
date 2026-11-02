# Shared ggplot2 typography and export dimensions.
# Export at the intended width; increasing raster DPI does not enlarge text.
infuse_plot_spec <- function() {
  list(width_mm = 170, height_mm = 145, family = "Liberation Sans",
       base_pt = 8, axis_title_pt = 8.5, tick_pt = 7.5,
       strip_pt = 8.5, tag_pt = 11, legend_pt = 7.5,
       value_pt = 6, key_mm = 2.3, png_dpi = 600)
}

theme_infuse_plot <- function(spec = infuse_plot_spec()) {
  ggplot2::theme_classic(base_size = spec$base_pt, base_family = spec$family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(size = spec$tick_pt, colour = "black"),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
      axis.title = ggplot2::element_text(size = spec$axis_title_pt, face = "bold"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = grid::unit(0.7, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = spec$strip_pt, face = "bold",
        margin = ggplot2::margin(t = 3, b = 3), lineheight = 0.95),
      panel.spacing = grid::unit(2.5, "mm"),
      plot.title = ggplot2::element_text(size = spec$strip_pt, face = "bold"),
      plot.tag = ggplot2::element_text(size = spec$tag_pt, face = "bold"),
      plot.margin = ggplot2::margin(3, 3, 3, 3),
      legend.position = "bottom", legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = spec$legend_pt),
      legend.title = ggplot2::element_text(size = spec$legend_pt, face = "bold"),
      legend.key.width = grid::unit(spec$key_mm, "mm"),
      legend.key.height = grid::unit(spec$key_mm, "mm"),
      legend.spacing.x = grid::unit(2.5, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(2, 0, 0, 0)
    )
}

# geom_text uses millimetres rather than theme element_text's points.
infuse_text_mm <- function(points) points / (72.27 / 25.4)
