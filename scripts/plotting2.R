library(tidyverse)
library(ggplot2)

# ── Load data ─────────────────────────────────────────────────────────────────
df <- read_csv("haplogroup_proportions_single_source_from_Diagnostics10_files.csv")

# ── Shared palette & theme tweaks ─────────────────────────────────────────────
damage_pal  <- c("none" = "#4E79A7", "low" = "#E8A838", "high" = "#C0392B")
seqtype_pal <- c("capture" = "#4E79A7", "shotgun" = "#59A14F")

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, colour = "grey40"),
    plot.caption       = element_text(size = 8, colour = "grey50")
  )

# ── Pre-processing ────────────────────────────────────────────────────────────
df <- df %>%
  mutate(
    simulated_reads = factor(simulated_reads, levels = sort(unique(simulated_reads))),
    damage          = factor(damage, levels = c("none", "low", "high")),
    correct         = true_haplogroup == predicted_haplogroup
  )

# Summarise accuracy (sum of proportion where prediction is correct)
acc <- df %>%
  group_by(species, seq_type, simulated_reads, damage) %>%
  summarise(pct_correct = sum(proportion[correct], na.rm = TRUE), .groups = "drop")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 – Accuracy by read depth & damage level, faceted by species × seq type
# ══════════════════════════════════════════════════════════════════════════════
p1 <- ggplot(acc,
             aes(x = simulated_reads, y = pct_correct,
                 colour = damage, group = damage)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.4) +
  facet_grid(species ~ seq_type) +
  scale_colour_manual(values = damage_pal, name = "Damage level") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title    = "Haplogroup prediction accuracy across read depths",
    subtitle = "Accuracy = summed proportion where predicted == true haplogroup",
    x        = "Simulated fragments",
    y        = "Proportion correctly predicted",
    caption  = "Rows = species  |  Columns = sequencing type"
  ) +
  base_theme

print(p1)
ggsave("plot1_accuracy_by_depth_damage.pdf", p1, width = 8, height = 5 * n_distinct(df$species))


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 – Capture vs shotgun accuracy, all read depths pooled
#           Boxplot + jitter coloured by damage
# ══════════════════════════════════════════════════════════════════════════════
p2 <- ggplot(acc,
             aes(x = seq_type, y = pct_correct, fill = seq_type)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = damage), width = 0.18, size = 2, alpha = 0.8) +
  facet_wrap(~ species, nrow = 1) +
  scale_fill_manual(values   = seqtype_pal, guide = "none") +
  scale_colour_manual(values = damage_pal, name = "Damage level") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title    = "Capture vs shotgun: distribution of prediction accuracy",
    subtitle = "Each point = one read-depth × damage combination",
    x        = "Sequencing type",
    y        = "Proportion correctly predicted"
  ) +
  base_theme

print(p2)
ggsave("plot2_capture_vs_shotgun.pdf", p2, width = 4 * n_distinct(df$species), height = 5)


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3 – Heatmap of mean proportion assigned to each predicted haplogroup
#           Rows = predicted haplogroup, Cols = read depth
#           Facets = species × damage × seq_type
# ══════════════════════════════════════════════════════════════════════════════
heatmap_df <- df %>%
  group_by(species, seq_type, damage, simulated_reads, true_haplogroup, predicted_haplogroup) %>%
  summarise(prop = mean(proportion, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(heatmap_df,
             aes(x = simulated_reads, y = predicted_haplogroup, fill = prop)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_grid(
    species + seq_type ~ damage,
    labeller = labeller(damage = label_both)
  ) +
  scale_fill_gradient2(
    low      = "#F5F5F5",
    mid      = "#76B7B2",
    high     = "#1A5276",
    midpoint = 0.5,
    limits   = c(0, 1),
    labels   = scales::percent_format(accuracy = 1),
    name     = "Mean\nproportion"
  ) +
  labs(
    title    = "Haplogroup assignment heatmap",
    subtitle = "Diagonal = correct assignment; off-diagonal = misassignment",
    x        = "Simulated fragments",
    y        = "Predicted haplogroup"
  ) +
  base_theme +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid       = element_blank(),
    legend.position  = "right"
  )

print(p3)
ggsave("plot3_heatmap_assignment.pdf", p3,
       width  = 10,
       height = 4 * n_distinct(df$species) * n_distinct(df$seq_type))


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 – Effect of damage at each read depth (bar chart of accuracy)
#           Good for a direct low/high/none comparison per condition
# ══════════════════════════════════════════════════════════════════════════════
p4 <- ggplot(acc,
             aes(x = damage, y = pct_correct, fill = damage)) +
  geom_col(width = 0.65, alpha = 0.85) +
  facet_grid(species ~ simulated_reads + seq_type,
             labeller = labeller(simulated_reads = label_both)) +
  scale_fill_manual(values = damage_pal, name = "Damage level") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title    = "Effect of damage on prediction accuracy",
    subtitle = "Faceted by fragment depth and sequencing type",
    x        = "Damage level",
    y        = "Proportion correctly predicted"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p4)
ggsave("plot4_damage_effect.pdf", p4,
       width  = 3 * n_distinct(df$simulated_reads),
       height = 4 * n_distinct(df$species))

