library(ggplot2)
library(ggtree)
library(dplyr)
library(ape)
library(ggrepel)
library(patchwork)
setwd("/Users/bfj994/Documents/PhD_related/costumGraphs/")
plot_mcmc_placements <- function(tree_file, input_file, title = NULL,
                                 y_buffer = 4, x_buffer_l = 0.002, x_buffer_r = 0.003) {
  
  # ------------------------------------------------------------------
  # 1. Read data
  # ------------------------------------------------------------------
  tree <- read.tree(tree_file)
  df   <- read.csv(input_file, header = TRUE, sep = "\t")
  colnames(df) <- make.names(colnames(df))
  df <- df %>%
    rename(
      source_node     = Source_1,
      log_likelihood  = Log.likelihood,
      branch_pos      = branch_position_derived,
      true_source     = set_label,
      branch_distance = branch_distance
    )
  
  # ------------------------------------------------------------------
  # 2. Frequency per placement position
  # ------------------------------------------------------------------
  df <- df %>%
    group_by(source_node, branch_pos) %>%
    mutate(freq = n()) %>%
    ungroup()
  
  # ------------------------------------------------------------------
  # 3. ggtree layout
  # ------------------------------------------------------------------
  p0     <- ggtree(tree)
  d_tree <- ggplot_build(p0)$data[[1]]
  nameList <- c(tree$tip.label, tree$node.label)
  
  # ------------------------------------------------------------------
  # 4. Map MCMC rows to branch coordinates
  # ------------------------------------------------------------------
  matched_idx <- match(df$source_node, nameList)
  df$plot_x <- NA_real_
  df$plot_y <- NA_real_
  valid <- !is.na(matched_idx)
  df$plot_x[valid] <- d_tree[matched_idx[valid], "x"] +
    df$branch_pos[valid] *
    (d_tree[matched_idx[valid], "xend"] - d_tree[matched_idx[valid], "x"])
  df$plot_y[valid] <- d_tree[matched_idx[valid], "y"]
  
  # Y jitter by log-likelihood
  median_ll <- median(df$log_likelihood, na.rm = TRUE)
  max_diff  <- max(abs(df$log_likelihood - median_ll), na.rm = TRUE)
  if (max_diff > 0) {
    df$plot_y <- df$plot_y +
      ifelse(df$log_likelihood < median_ll, -1, 1) *
      (abs(df$log_likelihood - median_ll) / max_diff) * 0.2
  }
  
  df_u <- df %>%
    select(source_node, plot_x, plot_y, log_likelihood, freq, true_source) %>%
    distinct() %>%
    mutate(log_freq = log1p(freq))
  
  # ------------------------------------------------------------------
  # 5. True source tip marker
  # ------------------------------------------------------------------
  true_sources    <- unique(df_u$true_source)
  true_source_idx <- match(true_sources, nameList)
  true_source_pts <- data.frame(
    x     = d_tree[true_source_idx, "xend"],
    y     = d_tree[true_source_idx, "y"],
    label = true_sources
  )
  
  # ------------------------------------------------------------------
  # 6. Zoom window bounds
  # ------------------------------------------------------------------
  xlim_zoom <- c(
    min(c(df_u$plot_x, true_source_pts$x), na.rm = TRUE) - x_buffer_l,
    max(c(df_u$plot_x, true_source_pts$x), na.rm = TRUE) + x_buffer_r
  )
  ylim_zoom <- c(
    min(c(df_u$plot_y, true_source_pts$y), na.rm = TRUE) - y_buffer,
    max(c(df_u$plot_y, true_source_pts$y), na.rm = TRUE) + y_buffer
  )
  
  # ------------------------------------------------------------------
  # Shared point layers
  # ------------------------------------------------------------------
  placement_pts <- geom_point(
    data = df_u,
    aes(x = plot_x, y = plot_y, color = log_likelihood, size = log_freq),
    alpha = 0.6, inherit.aes = FALSE
  )
  true_src_tri <- geom_point(
    data = true_source_pts,
    aes(x = x, y = y),
    shape = 17, size = 4, color = "black", inherit.aes = FALSE
  )
  shared_scales <- list(
    scale_color_gradientn(
      colors = c("grey70", "yellow", "orange", "red"),
      name   = "Log-likelihood"
    ),
    scale_size_continuous(
      name  = "Frequency [log1p]",
      range = c(1, 5)
    ),
    guides(alpha = "none")
  )
  
  # ------------------------------------------------------------------
  # 7. Full tree panel — xlim set to just past the rectangle's right edge
  # ------------------------------------------------------------------
  x_max_full <- max(d_tree$xend, na.rm = TRUE)
  # Add a fixed padding so the rect is never clipped
  x_right_full <- max(xlim_zoom[2], x_max_full) + 0.0005
  
  p_full <- p0 +
    coord_cartesian(xlim = c(0, x_right_full), clip = "off") +
    annotate("rect",
             xmin = xlim_zoom[1], xmax = xlim_zoom[2],
             ymin = ylim_zoom[1], ymax = ylim_zoom[2],
             fill = NA, color = "#E63946", linewidth = 0.8, linetype = "dashed"
    ) +
    placement_pts +
    true_src_tri +
    shared_scales +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "none",
      panel.background = element_blank(),
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      axis.line.x      = element_line(color = "black"),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank()
    )
  
  # ------------------------------------------------------------------
  # 8. Zoomed panel — coord_cartesian extended on the right for stretch
  # ------------------------------------------------------------------
  x_right_zoom <- xlim_zoom[2] + 0.001   # increase this for more stretch
  
  p_zoom <- p0 +
    geom_tippoint(size = 1.5, color = "grey40") +
    placement_pts +
    true_src_tri +
    geom_label_repel(
      data = true_source_pts,
      aes(x = x, y = y, label = label),
      size = 5, nudge_x = 0.001, box.padding = 0.3,
      segment.color = "grey50", inherit.aes = FALSE
    ) +
    shared_scales +
    coord_cartesian(
      xlim = c(xlim_zoom[1], x_right_zoom),
      ylim = ylim_zoom,
      clip = "on"
    ) +
    theme(
      legend.position = "none",
      panel.border    = element_blank()
    )
  
  # ------------------------------------------------------------------
  # 9. Combine with patchwork
  # ------------------------------------------------------------------
  combined <- (p_full + p_zoom) +
    plot_layout(widths = c(1, 1), guides = "collect")
  return(combined)
}

# ------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------
p <- plot_mcmc_placements(
  tree_file  = "Beluga.new.dnd",
  input_file = "BelugaSimRes/Beluga_Capture_1000x_5LowResult11_with_distances.tsv"
)
p
ggsave("beluga_mcmc_placements_1000.png", p, width = 12, height = 8, dpi = 600)

p5 <- plot_mcmc_placements(
  tree_file  = "Beluga.new.dnd",
  input_file = "BelugaSimRes/Beluga_Capture_50x_5LowResult11_with_distances.tsv"
)
p5
ggsave("beluga_mcmc_placements_50.png", p5, width = 12, height = 8, dpi = 600)

final_plot <- (p / p5) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )
ggsave(
  "combined_mcmc_panels.png",
  final_plot,
  width = 10,
  height = 13,
  dpi = 600
)

