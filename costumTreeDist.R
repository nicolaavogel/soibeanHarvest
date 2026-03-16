setwd("/Users/bfj994/Documents/PhD_related/costumGraphs/")

library(ape)
library(ggtree)
library(tidyverse)

# ---- Distance computation function ----
compute_label_distance <- function(label, set_label, layout, pos, tree) {
  label_nodes <- layout$node[layout$label == label]
  set_nodes   <- layout$node[layout$label == vg ]
  
  if (length(label_nodes) == 0 || length(set_nodes) == 0) {
    warning("Label or set_label not found in layout.")
    return(NA)
  }
  if (length(label_nodes) > 1 || length(set_nodes) > 1) {
    warning("Multiple matches found for label or set_label. Using first match.")
  }
  
  label_node <- label_nodes[1]
  set_node <- set_nodes[1]
  
  if (is.na(pos)) {
    warning("Position (pos) is NA.")
    return(NA)
  }
  
  if (label_node == set_node) {
    branch_length <- layout$branch.length[layout$node == label_node][1]
    label_pos <- branch_length * pos
    set_pos   <- branch_length
    dist <- abs(label_pos - set_pos)
    cat(sprintf("     - Same branch\n     - Branch length: %.4f | Label pos: %.4f | Set pos: %.4f | Distance: %.4f\n",
                branch_length, label_pos, set_pos, dist))
    return(dist)
  }
  
  path_nodes <- ape::nodepath(tree, from = label_node, to = set_node)
  branch_lengths <- sapply(path_nodes[-1], function(node_id) {
    layout$branch.length[layout$node == node_id][1]
  })
  
  first_branch_node <- path_nodes[2]
  first_branch_length <- layout$branch.length[layout$node == first_branch_node][1]
  partial_start <- first_branch_length * (1 - pos)
  
  total_dist <- sum(branch_lengths, na.rm = TRUE) - first_branch_length + partial_start
  
  cat(sprintf("     - Different branches\n     - Path nodes: %s\n     - Branch lengths: %s\n     - Partial start: %.4f | Total distance: %.4f\n",
              paste(path_nodes, collapse = " → "),
              paste(round(branch_lengths, 4), collapse = ", "),
              partial_start, total_dist))
  return(total_dist)
}

# ---- Tree + layout ----
tree_file <- "tree_dir/harpSeal.new.dnd"
tree <- read.tree(tree_file)
layout <- ggtree(tree)$data
set_label <- "KP942574.1_Phoca_groenlandica_isolate_Pag015_W691_A4_mitochondrion_partial_genome"

# ---- Batch process all matching .mcmc files ----
input_files <- Sys.glob("harpSeal/harpSeal1*Result1*.mcmc")

for (input_file in input_files) {
  cat("=== Processing:", input_file, "===\n")
  red_data <- read_tsv(input_file, show_col_types = FALSE)
  
  if (!all(c("Source_1", "branch_position_derived") %in% names(red_data))) {
    warning("Missing required columns in", input_file)
    next
  }
  
  red_data <- red_data %>%
    mutate(
      set_label = Source_1,
      branch_distance = map2_dbl(Source_1, branch_position_derived,
                                 ~ compute_label_distance(.x, set_label, layout, .y, tree))
    )
  
  output_file <- sub("\\.mcmc$", "_with_distances.tsv", input_file)
  write_tsv(red_data, output_file)
  cat("✓ Saved:", output_file, "\n\n")
}

############### PART 2 #######################
##### READING IN TSV FILES FOR PLOTTING ######


# ---- Parse metadata from filenames ----
input_files <- Sys.glob("BelugaSimRes//*_with_distance.tsv")

# ---- Summarize mean & median branch distances ----
summary_table <- map_dfr(input_files, function(file) {
  data <- read_tsv(file, show_col_types = FALSE)
  
  mean_dist   <- mean(data$branch_distance, na.rm = TRUE)
  median_dist <- median(data$branch_distance, na.rm = TRUE)
  
  tibble(
    file = basename(file),
    mean_branch_distance = mean_dist,
    median_branch_distance = median_dist
  )
})

################################################################################
### BELUGA SUMMARY TSV 

beluga <- read.csv("BelugaRes/summary_branch_distances.tsv", sep = '\t')

summary_table2 <- beluga %>%
  mutate(
    filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
    
    parsed = str_match(
      filename_no_ext,
      "^([A-Za-z]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
    ),
    
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

# Fix the mean_branch_distance column
df <- summary_table2 %>%
  mutate(mean_branch_distance = as.numeric(mean_branch_distance))|>
  filter(!is.na(method))

method_labels <- c(
  Capture = "Simulated Capture data (avg. fragment length 85bp)",
  Shotgun = "Simulated Shotgun data (avg. fragment length 45bp)"
)

#Summarize data: mean and 95% CI for each group
summary_df <- df %>%
  group_by(method, reads, simulation) %>%
  summarise(
    mean_mbd = mean(mean_branch_distance, na.rm = TRUE),
    sd_mbd = sd(mean_branch_distance, na.rm = TRUE),
    n = n(),
    se = sd_mbd / sqrt(n),
    ci95 = se * 1.96  # 95% CI assuming normality
  ) %>%
  ungroup()

# Plot: mean with error bars
beluga_plot <- ggplot(summary_df, aes(x = factor(reads), y = mean_mbd, color = simulation, group = simulation)) +
  geom_point(position = position_dodge(0.3), size = 4) +
  geom_errorbar(aes(ymin = mean_mbd - ci95, ymax = mean_mbd + ci95), 
                width = 0.5, position = position_dodge(0.3)) +
  labs(
    x = "Number of simulated fragments",
    y = "Mean Branch Distance with 95% confidence intervals",
    color = "Simulated Damage",
    title = "Beluga simulations"
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~method, labeller = labeller(method = method_labels))
beluga_plot

################################################################################
## harp seal
harp <- read.csv("harpSealRes/summary_branch_distances.tsv", sep = '\t')


summary_table2 <- harp %>%
  mutate(
    filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
    
    parsed = str_match(
      filename_no_ext,
      "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
      
    ),
    
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

df <- summary_table2 %>%
  mutate(mean_branch_distance = as.numeric(mean_branch_distance)) %>%
  filter(!is.na(method))

method_labels <- c(
  Capture = "Capture (avg. fragment length 85bp)",
  Shotgun = "Shotgun (avg. fragment length 45bp)"
)

harp_plot <- ggplot(df, aes(x = factor(reads), y = median_branch_distance, fill = simulation)) +
  geom_boxplot() +
  labs(
    x = "Number of simulated fragments",
    y = "Mean Branch Distance",
    fill = "Simulated Damage",
    title = "Harp Seal simulations"
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~method, labeller = labeller(method = method_labels))
harp_plot

################################################################################
## Narwhal
nar <- read.csv("NarwhalRes/summary_branch_distances.tsv", sep = '\t')


summary_table2 <- nar %>%
  mutate(
    filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
    
    parsed = str_match(
      filename_no_ext,
      "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
      
    ),
    
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

df <- summary_table2 %>%
  mutate(mean_branch_distance = as.numeric(mean_branch_distance)) %>%
  filter(!is.na(method))

method_labels <- c(
  Capture = "Capture (avg. fragment length 85bp)",
  Shotgun = "Shotgun (avg. fragment length 45bp)"
)

nar_plot <- ggplot(df, aes(x = factor(reads), y = median_branch_distance, fill = simulation)) +
  geom_boxplot() +
  labs(
    x = "Number of simulated fragments",
    y = "Mean Branch Distance",
    fill = "Simulated Damage",
    title = "Narwhal simulations"
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~method, labeller = labeller(method = method_labels))
nar_plot

################################################################################
## harp seal
bow <- read.csv("bowHeadWhaleRes//summary_branch_distances.tsv", sep = '\t')


summary_table2 <- bow %>%
  mutate(
    filename_no_ext = str_remove(file, "_with_distances\\.tsv$"),
    
    parsed = str_match(
      filename_no_ext,
      "^([A-Za-z0-9]+)_([A-Za-z]+)_([0-9]+)x_([0-9]+)([A-Za-z]+)Result([0-9]+)$"
      
    ),
    
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

df <- summary_table2 %>%
  mutate(mean_branch_distance = as.numeric(mean_branch_distance)) %>%
  filter(!is.na(method))

method_labels <- c(
  Capture = "Capture (avg. fragment length 85bp)",
  Shotgun = "Shotgun (avg. fragment length 45bp)"
)

bow_plot <- ggplot(df, aes(x = factor(reads), y = median_branch_distance, fill = simulation)) +
  geom_boxplot() +
  labs(
    x = "Number of simulated fragments",
    y = "Mean Branch Distance",
    fill = "Simulated Damage",
    title = "Bowhead Whale simulations"
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~method, labeller = labeller(method = method_labels))
bow_plot
