setwd("/Users/bfj994/Documents/PhD_related/costumGraphs/")
library(ggplot2)
library(tidyverse)
library(stringr)
library(ape)
library(ggtree)
library(dplyr)
library(readxl)


plot_simulation_summary <- function(species_name, filepath, method_labels = NULL, extended_species = FALSE) {
  # Choose regex based on flag
  regex_pattern <- if (extended_species) {
    "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  } else {
    "^([A-Za-z]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  }
  
  # Read input
  df_raw <- read.csv(filepath, sep = '\t')
  
  summary_table <- df_raw %>%
    mutate(
      filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
      
      parsed = str_match(filename_no_ext, regex_pattern),
      
      species     = parsed[,2],
      method      = parsed[,3],
      reads       = as.numeric(parsed[,4]),
      bootstrap   = as.integer(parsed[,5]),
      simulation  = parsed[,6],
      chain       = parsed[,7],
      coverage    = round(
        if_else(method == "Capture",
                (reads * 85) / 16500,
                (reads * 45) / 16500
        ), 2
      ),
      chain_label = paste0("1.", chain)
    ) %>%
    select(-parsed, -filename_no_ext)
  
  df <- summary_table %>%
    mutate(mean_branch_distance = as.numeric(mean_branch_distance)) %>%
    filter(!is.na(method))
  
  # Default method labels
  if (is.null(method_labels)) {
    method_labels <- c(
      Capture = "Simulated Capture (85bp)",
      Shotgun = "Simulated Shotgun (45bp)"
    )
  }
  
  # Summary
  summary_df <- df %>%
    group_by(method, reads, simulation) %>%
    summarise(
      mean_mbd = mean(mean_branch_distance, na.rm = TRUE),
      sd_mbd = sd(mean_branch_distance, na.rm = TRUE),
      n = n(),
      se = sd_mbd / sqrt(n),
      ci95 = se * 1.96,
      .groups = "drop"
    )
  
  # Plot
  p <- ggplot(summary_df, aes(x = factor(reads), y = mean_mbd, color = simulation, group = simulation)) +
    geom_point(position = position_dodge(0.3), size = 4) +
    geom_errorbar(aes(ymin = mean_mbd - ci95, ymax = mean_mbd + ci95), 
                  width = 0.5, position = position_dodge(0.3)) +
    labs(
      x = "Number of simulated fragments",
      y = "Mean Branch Distance with 95% CI",
      color = "Simulated Damage",
      title = paste(species_name, "Simulations")
    ) +
    theme_bw(base_size = 16) +
    facet_wrap(~method, labeller = labeller(method = method_labels), ncol = 1)
  
  return(p)
}


beluga_plot <- plot_simulation_summary(
  species_name = "Beluga",
  filepath = "BelugaRes/summary_branch_distances.tsv"
)
beluga_plot
ggsave("beluga_sim.png", beluga_plot, width = 12, height = 10)
harp_plot <- plot_simulation_summary(
  species_name = "Harp Seal",
  filepath = "harpSealRes/summary_branch_distances.tsv", extended_species = TRUE
)
harp_plot
ggsave("harp_sim.png", harp_plot, width = 12, height = 10)

nar_plot <- plot_simulation_summary(
  species_name = "Narwhal",
  filepath = "NarwhalRes/summary_branch_distances.tsv", extended_species = TRUE
)
nar_plot
ggsave("nar_sim.png", nar_plot, width = 12, height = 10)

bow_plot <- plot_simulation_summary(
  species_name = "Bowhead Whale",
  filepath = "bowHeadWhaleRes/summary_branch_distances.tsv", extended_species = TRUE
)
bow_plot
ggsave("bow_sim.png", bow_plot, width = 12, height = 10)

###############################################################################
#ALL IN ONE PLOTS 
library(tidyverse)
library(stringr)
library(patchwork)

plot_simulation_summary <- function(
    species_name, 
    filepath, 
    method_labels = NULL, 
    extended_species = FALSE,
    y_limits = NULL
) {
  
  # Select regex pattern
  regex_pattern <- if (extended_species) {
    "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  } else {
    "^([A-Za-z]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  }
  
  # Read input
  df_raw <- read.csv(filepath, sep = '\t')
  
  summary_table <- df_raw %>%
    mutate(
      filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
      parsed = str_match(filename_no_ext, regex_pattern),
      species     = parsed[,2],
      method      = parsed[,3],
      reads       = as.numeric(parsed[,4]),
      bootstrap   = as.integer(parsed[,5]),
      simulation  = parsed[,6],
      chain       = parsed[,7],
      coverage    = round(
        if_else(method == "Capture",
                (reads * 85) / 16500,
                (reads * 45) / 16500
        ), 2
      ),
      chain_label = paste0("1.", chain)
    ) %>%
    select(-parsed, -filename_no_ext)
  
  df <- summary_table %>%
    mutate(mean_branch_distance = as.numeric(mean_branch_distance)) %>%
    filter(!is.na(method))
  
  # Default method labels
  if (is.null(method_labels)) {
    method_labels <- c(
      Capture = "Simulated Capture (~85bp)",
      Shotgun = "Simulated Shotgun (~45bp)"
    )
  }
  
  # Summary stats
  summary_df <- df %>%
    group_by(method, reads, simulation) %>%
    summarise(
      mean_mbd = mean(mean_branch_distance, na.rm = TRUE),
      sd_mbd = sd(mean_branch_distance, na.rm = TRUE),
      n = n(),
      se = sd_mbd / sqrt(n),
      ci95 = se * 1.96,
      .groups = "drop"
    )
  
  # Plot
  p <- ggplot(summary_df, aes(x = factor(reads), y = mean_mbd, 
                              color = simulation, group = simulation)) +
    geom_point(position = position_dodge(0.3), size = 4) +
    geom_errorbar(aes(ymin = mean_mbd - ci95, ymax = mean_mbd + ci95),
                  width = 0.5, position = position_dodge(0.3)) +
    labs(
      x = "Number of simulated fragments",
      y = "Mean Branch Distance with 95% CI",
      color = "Simulated Damage",
      title = paste(species_name, "Simulations")
    ) +
    scale_color_manual(
      values = c(
        "High" = "#C62828",
        "Low"  = "#2E8B57",
        "None" = "#005F9E"
      )
    )+
    theme_bw(base_size = 14) +
    facet_wrap(~method, labeller = labeller(method = method_labels), ncol = 1)
  
  # Shared y-axis limits if provided
  if (!is.null(y_limits)) {
    p <- p + ylim(y_limits)
  }
  
  return(p)
}

# -----------------------------------------
# PART 2: Build combined plot with shared Y-scale
# -----------------------------------------
library(tidyverse)
library(stringr)
library(ggplot2)

# Filepaths
files <- list(
  Beluga       = "BelugaRes/summary_branch_distances.tsv",
  HarpSeal     = "harpSealRes/summary_branch_distances.tsv",
  Narwhal      = "NarwhalRes/summary_branch_distances.tsv",
  BowheadWhale = "bowHeadWhaleRes/summary_branch_distances.tsv"
)

# Function to read and parse one species
read_simulation_data <- function(species_name, filepath, extended_species = FALSE) {
  
  regex_pattern <- if (extended_species) {
    "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  } else {
    "^([A-Za-z]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  }
  
  df_raw <- read.csv(filepath, sep = "\t")
  
  df <- df_raw %>%
    mutate(
      filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
      parsed = str_match(filename_no_ext, regex_pattern),
      species_file = parsed[,2],
      method = parsed[,3],
      reads = as.numeric(parsed[,4]),
      bootstrap = as.integer(parsed[,5]),
      simulation = parsed[,6],
      chain = parsed[,7],
      species = species_name,
      mean_branch_distance = as.numeric(mean_branch_distance)
    ) %>%
    filter(!is.na(method)) %>%
    select(species, method, reads, bootstrap, simulation, chain, mean_branch_distance)
  
  return(df)
}

# Read all species into one dataframe
all_df <- bind_rows(
  read_simulation_data("Beluga", files$Beluga, extended_species = FALSE),
  read_simulation_data("Harp Seal", files$HarpSeal, extended_species = TRUE),
  read_simulation_data("Narwhal", files$Narwhal, extended_species = TRUE),
  read_simulation_data("Bowhead Whale", files$BowheadWhale, extended_species = TRUE)
)

# Optional: make labels prettier if needed
# Adjust these if your real labels differ
all_df <- all_df %>%
  mutate(
    simulation = recode(simulation,
                        "High" = "High",
                        "Low" = "Low",
                        "None" = "None"
    )
  )

method_labels <- c(
  Capture = "Simulated Capture (~85 bp)",
  Shotgun = "Simulated Shotgun (~45 bp)"
)

# Summarise
summary_df <- all_df %>%
  group_by(species, method, reads, simulation) %>%
  summarise(
    mean_mbd = mean(mean_branch_distance, na.rm = TRUE),
    sd_mbd = sd(mean_branch_distance, na.rm = TRUE),
    n = n(),
    se = sd_mbd / sqrt(n),
    ci95 = se * 1.96,
    .groups = "drop"
  )

# Global y-axis range including CI
global_y <- range(
  c(summary_df$mean_mbd - summary_df$ci95,
    summary_df$mean_mbd + summary_df$ci95),
  na.rm = TRUE
)

# Combined plot
combined_plot <- ggplot(
  summary_df,
  aes(x = factor(reads), y = mean_mbd, color = simulation, group = simulation)
) +
  geom_point(position = position_dodge(0.3), size = 4) +
  geom_errorbar(
    aes(ymin = mean_mbd - ci95, ymax = mean_mbd + ci95),
    width = 0.5,
    position = position_dodge(0.3)
  ) +
  labs(
    x = "Number of simulated fragments",
    y = "Mean Branch Distance with 95% CI",
    color = "Simulated Damage"
  ) +
  scale_color_manual(
    values = c(
      "High" = "#C62828",
      "Low"  = "#2E8B57",
      "None" = "#005F9E"
    )
  ) +
  facet_grid(species ~ method, labeller = labeller(method = method_labels)) +
  coord_cartesian(ylim = global_y) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

combined_plot

ggsave("all_species_sim.png", combined_plot, width = 16, height = 12, dpi = 600)
##############################################################################
## FULLY ONE 

read_simulation_file <- function(filepath, extended_species = FALSE) {
  
  regex_pattern <- if (extended_species) {
    "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  } else {
    "^([A-Za-z]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
  }
  
  df_raw <- read.csv(filepath, sep = '\t')
  
  df <- df_raw %>%
    mutate(
      filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
      parsed = str_match(filename_no_ext, regex_pattern),
      
      species     = parsed[,2],
      method      = parsed[,3],
      reads       = as.numeric(parsed[,4]),
      bootstrap   = as.integer(parsed[,5]),
      simulation  = parsed[,6],
      chain       = parsed[,7],
      
      coverage    = round(
        if_else(method == "Capture",
                (reads * 85) / 16500,
                (reads * 45) / 16500),
        2
      ),
      
      mean_branch_distance = as.numeric(mean_branch_distance)
    ) %>%
    select(-parsed, -filename_no_ext) %>%
    filter(!is.na(method))
  
  return(df)
}
files <- list(
  Beluga       = list(path = "BelugaRes/summary_branch_distances.tsv", extended = FALSE),
  HarpSeals    = list(path = "harpSealRes/summary_branch_distances.tsv", extended = TRUE),
  Narwhal      = list(path = "NarwhalRes/summary_branch_distances.tsv", extended = TRUE),
  BowheadWhale = list(path = "bowHeadWhaleRes/summary_branch_distances.tsv", extended = TRUE)
)

all_df <- purrr::map_dfr(names(files), function(sp) {
  info <- files[[sp]]
  read_simulation_file(info$path, info$extended) %>%
    mutate(species_name = sp)
})

summary_df <- all_df %>%
  group_by(species_name, method, reads, simulation) %>%
  summarise(
    mean_mbd = mean(mean_branch_distance, na.rm = TRUE),
    sd_mbd = sd(mean_branch_distance, na.rm = TRUE),
    n = n(),
    se = sd_mbd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  )
library(ggplot2)

p <- ggplot(summary_df,
            aes(x = factor(reads), y = mean_mbd,
                color = simulation, group = simulation)) +
  
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_errorbar(aes(ymin = mean_mbd - ci95,
                    ymax = mean_mbd + ci95),
                width = 0.3, position = position_dodge(0.3)) +
  
  facet_grid(species_name ~ method,
             labeller = labeller(method = c(
               Capture = "Simulated Capture (~85bp)",
               Shotgun = "Simulated Shotgun (~45bp)"
             )),
             scales = "free_y") +
  scale_color_manual(
    values = c(
      "High" = "#C62828",
      "Low"  = "#2E8B57",
      "None" = "#005F9E"
    )
  )+
  labs(
    x = "Number of simulated fragments",
    y = "Mean Branch Distance ± 95% CI",
    color = "Simulated Damage"
  ) +
  
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14),
    legend.position = "bottom"
  )

p
ggsave("combination_single_source_sim_free_y.png", p, width = 14, height = 8)

###############################################################################
library(ape)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggtree)
library(tidyr)

# Clade colors used across figures
legend_colors <- c(
  "A" = "#005F9E",
  "B" = "#C98900",
  "C" = "#2E8B57",
  "D" = "#6B2F8A",
  "E" = "#C62828"
)

# Read tree and metadata
tree <- read.tree("harpSeal.new.dnd")

clade_info <- read_excel("ref_mitogenome_clade_info.xlsx") |>
  filter(species == "harp seal")

# Pairwise tip distances
dist_matrix <- cophenetic(tree)

# Keep only samples present in both sources
tips_in_tree <- intersect(rownames(dist_matrix), clade_info$fasta_header)

dist_matrix <- dist_matrix[tips_in_tree, tips_in_tree]
clade_info <- clade_info |>
  filter(fasta_header %in% tips_in_tree)

# Keep clades in a stable order
clade_list <- sort(unique(clade_info$clade))

# Compute average pairwise distances between clades
avg_dist_clade <- matrix(
  NA_real_,
  nrow = length(clade_list),
  ncol = length(clade_list),
  dimnames = list(clade_list, clade_list)
)

for (clade1 in clade_list) {
  for (clade2 in clade_list) {
    tips1 <- clade_info$fasta_header[clade_info$clade == clade1]
    tips2 <- clade_info$fasta_header[clade_info$clade == clade2]
    
    dists <- dist_matrix[tips1, tips2, drop = FALSE]
    avg_dist_clade[clade1, clade2] <- mean(dists, na.rm = TRUE)
  }
}

# ------------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------------

df_plot <- as.data.frame(as.table(avg_dist_clade))
colnames(df_plot) <- c("Clade1", "Clade2", "AvgDistance")
df_plot <- df_plot |>
  mutate(
    Clade1 = factor(Clade1, levels = clade_list),
    Clade2 = factor(Clade2, levels = rev(clade_list))
  )

heatmap_plot <- ggplot(df_plot, aes(x = Clade1, y = Clade2, fill = AvgDistance)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(AvgDistance, 4)), size = 4) +
  scale_fill_gradient(
    low = "white",
    high = "#005F9E",
    name = "Mean\ndistance"
  ) +
  labs(
    title = "Average Pairwise Distances Between Beluga Clades",
    x = "Clade",
    y = "Clade"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13, color = "black"),
    panel.grid = element_blank()
  )

heatmap_plot
ggsave("HarpSeal_clade_distance_heatmap.png", heatmap_plot, width = 7, height = 6, dpi = 600)

# ------------------------------------------------------------------
# Neighbor-joining tree of clades
# ------------------------------------------------------------------

clade_dist <- as.dist(avg_dist_clade)
clade_tree <- nj(clade_dist)

# Use only clade colors that are present in the clade tree
tip_colors <- legend_colors[clade_tree$tip.label]

# Build ggtree object once
p_tree <- ggtree(clade_tree, layout = "dendrogram", linewidth = 1)

# Extract plotting data and add tip colors
tree_data <- p_tree$data |>
  mutate(
    tip_col = ifelse(isTip, label, NA_character_)
  )

p0 <- ggtree(clade_tree)
x_max <- max(p0$data$x, na.rm = TRUE)

clade_tree_plot <- p0 %<+% tree_data +
  geom_tree(linewidth = 1, color = "grey70") +
  geom_tiplab(
    aes(color = label),
    hjust = -0.2,
    size = 5,
    fontface = "bold"
  ) +
  geom_tippoint(aes(color = label), size = 3.5) +
  geom_text2(
    aes(label = round(branch.length, 5)),
    size = 3.8,
    color = "black",
    hjust = -0.35,
    vjust = -0.5
  ) +
  scale_color_manual(values = tip_colors, guide = "none") +
  theme_tree2() +
  xlim(0, x_max +0.001) +
  labs(title = "Tree of the Narwhal Haplogroups") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

clade_tree_plot
ggsave("HarpSeal_hap_tree.png", clade_tree_plot, width = 12, height = 5, dpi = 600)

################################################################################
# Two source plots 

setwd("/Users/bfj994/Documents/PhD_related/costumGraphs/Beluga_2_res/")

tree <- read.tree("../Beluga.new.dnd")
true_source_R1 <- "DL0852"
true_sources_R2 <- c("DL0852", "DL1661")
clade_info <- read_excel("../ref_mitogenome_clade_info.xlsx") |>
  filter(species == "beluga")


B50R1 <- read.csv("B50xR1.txt", header = FALSE)
B50R2 <- read.csv("B50R2.txt", header = FALSE, sep = "\t")
B500R1 <- read.csv("B500xR1.txt", header = FALSE)
B500R2 <- read.csv("B500R2.txt", header = FALSE, sep = "\t")
B100R1 <- read.csv("B100xR1.txt", header = FALSE)
B100R2 <- read.csv("B100R2.txt", header = FALSE, sep = "\t")
B1000R1 <- read.csv("B1000xR1.txt", header = FALSE)
B1000R2 <- read.csv("B1000R2.txt", header = FALSE, sep = "\t")

library(ggtree)
library(dplyr)
library(patchwork)  # or cowplot

make_tree_plot <- function(tree, clade_info, true_sources, pred_sources, title = NULL) {
  # Build metadata
  all_nodes <- data.frame(
    node = 1:(Ntip(tree) + tree$Nnode),
    label = c(tree$tip.label, tree$node.label)
  )
  
  meta <- all_nodes %>%
    left_join(clade_info, by = c("label" = "fasta_header")) %>%
    mutate(
      source = case_when(
        label %in% true_sources ~ "True Source",
        label %in% pred_sources ~ "Predicted Source",
        TRUE ~ "NA_character_"
      )
    )
  
  meta$branch_color <- ifelse(meta$source == "Predicted Source", "Predicted Source", meta$clade)
  
  legend_colors <- c(
    "A" = "#005F9E",
    "B" = "#C98900",
    "C" = "#2E8B57",
    "D" = "#6B2F8A",
    "Predicted Source" = "#C62828"
  )
  legend_shapes <- c("True Source" = 17)
  
  # Plot
  p <- ggtree(tree) %<+% meta +
    geom_tree(aes(color = branch_color)) +
    geom_tippoint(aes(shape = source), size = 3, na.rm = TRUE, show.legend = TRUE) +
    scale_color_manual(
      values = legend_colors,
      na.value = "grey85",
      breaks = names(legend_colors),
      name = "Haplogroup"
    )+
    scale_shape_manual(values = legend_shapes, name = NULL, na.translate = FALSE) +
    guides(
      color = guide_legend(override.aes = list(linetype = 1, shape = NA)),
      shape = guide_legend(override.aes = list(color = "black"))
    ) +
      theme(
        text = element_text(size = 14),
        legend.position = "right"
      )+
    ggtitle(title)
  
  return(p)
}

p_B50R1   <- make_tree_plot(tree, clade_info, true_source_R1, as.character(B50R1$V1), title = "Beluga 100 reads - 1 source")
p_B500R1  <- make_tree_plot(tree, clade_info, true_source_R1, as.character(B500R1$V1), title = "Beluga 1000 reads - 1 source")
p_B100R1  <- make_tree_plot(tree, clade_info, true_source_R1, as.character(B100R1$V1), title = "Beluga 200 reads - 1 source")
p_B1000R1 <- make_tree_plot(tree, clade_info, true_source_R1, as.character(B1000R1$V1), title = "Beluga 2000 reads - 1 source")

p_B50R2   <- make_tree_plot(tree, clade_info, true_sources_R2, as.character(B50R2$V1), title = "Beluga 100 reads - 2 source")
p_B500R2  <- make_tree_plot(tree, clade_info, true_sources_R2, as.character(B500R2$V1), title = "Beluga 500 reads - 2 source")
p_B100R2  <- make_tree_plot(tree, clade_info, true_sources_R2, as.character(B100R2$V1), title = "Beluga 100 reads - 2 source")
p_B1000R2 <- make_tree_plot(tree, clade_info, true_sources_R2, as.character(B1000R2$V1), title = "Beluga 1000 reads - 2 source")

library(cowplot)
# Extract one shared legend from one plot
shared_legend <- get_legend(
  p_B50R1 + theme(legend.position = "bottom")
)

# Remove legends from individual plots
p_B50R1_noleg   <- p_B50R1   + theme(legend.position = "none")
p_B100R1_noleg  <- p_B100R1  + theme(legend.position = "none")
p_B500R1_noleg  <- p_B500R1  + theme(legend.position = "none")
p_B1000R1_noleg <- p_B1000R1 + theme(legend.position = "none")

p_B50R2_noleg   <- p_B50R2   + theme(legend.position = "none")
p_B100R2_noleg  <- p_B100R2  + theme(legend.position = "none")
p_B500R2_noleg  <- p_B500R2  + theme(legend.position = "none")
p_B1000R2_noleg <- p_B1000R2 + theme(legend.position = "none")

# Build row of plots, keeping titles
plot_R1_panels <- plot_grid(
  p_B50R1_noleg, p_B100R1_noleg, p_B500R1_noleg, p_B1000R1_noleg,
  nrow = 1,
  align = "h"
)

plot_R2_panels <- plot_grid(
  p_B50R2_noleg, p_B100R2_noleg, p_B500R2_noleg, p_B1000R2_noleg,
  nrow = 1,
  align = "h"
)

# Add shared legend at bottom
plot_R1 <- plot_grid(
  plot_R1_panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

plot_R2 <- plot_grid(
  plot_R2_panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)
ggsave("Beluga_2sources.png", plot = plot_R2, width = 22, heigh = 10, dpi = 600)
ggsave("Beluga_1sources.png", plot = plot_R1, width = 22, heigh = 10, dpi = 600)




