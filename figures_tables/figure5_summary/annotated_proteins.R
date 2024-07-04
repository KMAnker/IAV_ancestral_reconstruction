library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)

setwd("PATH_TO_PROJECT_DIRECTORY")

shap_positions <- read_table("figures_tables/figure4_S4_S5_xgboost/shap_positions.txt") %>%
  mutate(track = "shap")

meme_positions <- read.table("selection_analysis/meme_selection.txt", header = TRUE) %>%
  rename(segment = protein)  %>%
  mutate(track = "pos selection")

meme_positions <- meme_positions %>%
  filter(!segment == "PB1-F2")
meme_positions_PA_X <- meme_positions %>%
  filter(segment == "PA-X")

segment_lengths <- c("PB2" = 759, "PB1" = 757, "PB1-F2" = 101, "PA" = 716, "PA-X" = 260, 
                     "H1" = 565, "H3" = 566, "NP" = 498, "N1" = 469, "N2" = 469, "M1" = 252, 
                     "M2" = 97, "NS1" = 237, "NEP" = 121)

segment_lengths_df <- as.data.frame(segment_lengths) %>%
  rownames_to_column("segment") %>%
  rename(length = segment_lengths)

combined <- bind_rows(meme_positions, shap_positions)


plotting_df <- combined %>%
  left_join(segment_lengths_df, by = "segment")
plotting_df$segment <- factor(plotting_df$segment, levels=c("PB2", "PB1", "PB1-F2", "PA",
                                                                  "PA-X", "H1", "H3", "NP",
                                                                  "N1", "N2", "M1", "M2", 
                                                                  "NS1", "NEP"))
plotting_df$trait <- factor(plotting_df$trait, levels=c("human-human", "swine-swine", "human-swine", "swine-human"))

domains <- read_table("figures_tables/figure5_summary/protein_domains.txt")
domains$segment <- factor(domains$segment, levels=c("PB2", "PB1", "PB1-F2", "PA",
                                                            "PA-X", "H1", "H3", "NP",
                                                            "N1", "N2", "M1", "M2", 
                                                            "NS1", "NEP"))
calculate_midpoints <- function(df) {
  df %>%
    group_by(segment) %>%
    mutate(midpoint = case_when(
      row_number() == n() ~ (start + stop) / 2, # Last row of the segment
      stop >= lead(start, default = Inf) ~ (start + lead(start)) / 2, # stop >= start of next row
      stop < lead(start, default = Inf) ~ (start + stop) / 2          # stop < start of next row
    )) %>%
    ungroup()
}
domains_plotting_df <- calculate_midpoints(domains)
domains_order <- names(sort(table(domains_plotting_df$domain),decreasing = TRUE))
domains_plotting_df$domain <- factor(domains_plotting_df$domain, levels = domains_order)

overlapping_positions <- plotting_df %>%
  group_by(segment, position) %>%
  filter(n_distinct(track) > 1) %>%
  ungroup() %>%
  select(segment, position) %>%
  distinct(position, .keep_all = TRUE)

# Create a data frame for custom y-axis labels
y_labels_df <- data.frame(
  track = unique(plotting_df$track),
  y = unique(plotting_df$track)
)

# create dataframe for marking alternative reading frames
alt_orfs_marks <- data.frame(
  segment = c("PB1", "PB1", "PA", "PA", "PA-X", "M1", "M1", "M1", "M2", "NS1", "NS1", "NS1", "NEP"),
  position = c(31, 132, 1, 260, 190, 1, 9, 238, 9, 1, 10, 167, 10)                   
)
alt_orfs_marks$segment <- factor(alt_orfs_marks$segment, levels=c("PB2", "PB1", "PB1-F2", "PA",
                                                    "PA-X", "H1", "H3", "NP",
                                                    "N1", "N2", "M1", "M2", 
                                                    "NS1", "NEP"))
annotated_proteins <- ggplot(plotting_df) +
  # Draw the gene as a colored bar for each segment with functional domains
  geom_rect(aes(xmin = 0, xmax = length, ymin = 0, ymax = 0.8), fill = "lightgrey") +
  geom_rect(data = domains_plotting_df, aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.8, fill = domain)) +
  geom_text(data = domains_plotting_df, aes(x = midpoint, y = 0.4, label = domain), size = 2) +
  scale_fill_manual(values = c("grey60", "#a9b4a3", "#d9ae94", "#c0d4d8", "#d79b73", "#ccd5c3", 
                               "#cec8b4", "#e3ba9fcc", "#c0d4d8","#60aaa9bb", "#e8d1c9", "#ccd5a7ff", 
                               "#bfcfce","#d4b9a8", "#e4d4d1", "#d4c4b8ff", "#d8cfc1","#cfd5a1fa", 
                               "#bfa09c", "#92aba1ee", "#b7c3c8", "#d1e0df", "#ccd1c6", "#ccd5c3ff", 
                               "#d9d1c5", "#d6c3c3", "#ccd5d1dd","#e3ba9f", "#a9a9b8", "#dcd0c2",
                               "#dcbfbe", "#dcbeca","#bab2a8", "#ccd0c2ff"), guide = "none") +
  # Add points for the features
  geom_point(aes(x = position, y = track, color = trait, shape = track), size = 1.5, alpha = 0.85, stroke = 0.5) +
  guides(shape = "none") +
  scale_color_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  scale_shape_manual(values = c(16,17)) +
  #8,19,20,24,73
  # Add vertical lines for overlapping positions
  geom_vline(data = overlapping_positions, aes(xintercept = position), color = "grey40", linewidth = 0.1) +
  # Add text labels for the overlapping positions
  geom_text_repel(data = overlapping_positions, aes(x = position, y = 2.5, label = position),
            size = 1.8, color = "black", box.padding = 0.05, min.segment.length = 0.01, point.padding = 0.01, 
            vjust = 1.1, segment.color = NA) +
  # Add custom y-axis labels with shapes
  geom_point(data = y_labels_df, aes(x = -5, y = track, shape = track), color = "grey", size = 1, inherit.aes = FALSE) +
  geom_text(data = y_labels_df, aes(x = -7.8, y = track, label = track), size = 2, hjust = 1.04, inherit.aes = FALSE) +
  # Add marks for alternative ORFs
  geom_point(data = alt_orfs_marks, aes(x = position, y = 0), shape = 73, color = "grey40", size = 1, inherit.aes = FALSE) +
  # Customize the plot
  labs(title = "",
       x = "Amino acid position",
       y = "") +
  theme_void() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = unit(7, "points")),
    legend.key.size = unit(6, "points"),
    legend.title = element_blank(),
    axis.title.x = element_text(size = 6, margin = margin(t = 10)),
    axis.text.x = element_text(size = unit(5, "points"), vjust = -1),
    axis.ticks.x = element_line(size = 0.5, color = "darkgrey"),
    axis.ticks.length.x  = unit(3, "points"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = unit(7, "points"), hjust = -0.01),
    panel.spacing.y = unit(0.8, "lines"),
    panel.grid.major.x = element_line(color = "darkgrey",size = 0.05,linetype = 5),
    panel.grid.minor.x = element_line(color = "lightgrey",size = 0.05,linetype = 5),
    plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  scale_x_continuous(n.breaks = 8, minor_breaks = seq(0, 759, by=10), expand = expansion(mult = c(0.09, 0.04))) +
  scale_y_discrete(expand = c(-2.2, 0.3)) +
  # Facet the plot based on the segment
  facet_wrap(~ segment, scales = "free_x", ncol = 1, strip.position = "left")



ggsave("figures_tables/figure5_summary/annotated_proteins.png", plot = annotated_proteins, dpi = 300, width = 18 , height = 25, units = "cm", bg = "white")

