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
  dplyr::rename(segment = protein)  %>%
  mutate(track = "pos selection")

meme_positions <- meme_positions %>%
  filter(!segment == "PB1-F2")

bayes_positions <- read_table("figures_tables/figure2_mut_count_plots/bayes_positions.txt") %>%
  mutate(track = "mut rate")

segment_lengths <- c("PB2" = 759, "PB1" = 757, "PB1-F2" = 101, "PA" = 716, "PA-X" = 260, 
                     "H1" = 565, "H3" = 566, "NP" = 498, "N1" = 469, "N2" = 469, "M1" = 252, 
                     "M2" = 97, "NS1" = 237, "NEP" = 121)

segment_lengths_df <- as.data.frame(segment_lengths) %>%
  rownames_to_column("segment") %>%
  dplyr::rename(length = segment_lengths)

combined <- bind_rows(meme_positions, shap_positions, bayes_positions)


plotting_df <- combined %>%
  left_join(segment_lengths_df, by = "segment")
plotting_df$segment <- factor(plotting_df$segment, levels=c("PB2", "PB1", "PB1-F2", "PA",
                                                                  "PA-X", "H1", "H3", "NP",
                                                                  "N1", "N2", "M1", "M2", 
                                                                  "NS1", "NEP"))
plotting_df$trait <- factor(plotting_df$trait, levels=c("human-human", "swine-swine", "human-swine", "swine-human"))
plotting_df$track <- factor(plotting_df$track, levels=c("pos selection", "mut rate", "shap"))
                              
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
overlapping_positions$segment <- factor(overlapping_positions$segment, levels = c("PB2", "PB1", "PB1-F2", "PA", "PA-X",
                                                                                  "H1", "H3", "NP", "N1", "N2", "M1", "M2",
                                                                                  "NS1", "NEP"))


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
  geom_point(aes(x = position, y = track, color = trait, shape = track), size = 1.5, alpha = 0.75, stroke = 0.5) +
  guides(shape = "none") +
  scale_color_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  scale_shape_manual(values = c(19,15,17)) +
  # Add vertical lines for overlapping positions
  geom_vline(data = overlapping_positions, aes(xintercept = position), color = "grey40", linewidth = 0.1) +
  # Add text labels for the overlapping positions
  geom_text_repel(data = overlapping_positions, aes(x = position, y = 4, label = position),
            size = 1.7, color = "black", box.padding = 0.05, min.segment.length = 0.01, point.padding = 0.01, 
            vjust = 1, segment.color = NA) +
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
  scale_y_discrete(expand = c(-1.2, 0)) +
  # Facet the plot based on the segment
  facet_wrap(~ segment, scales = "free_x", ncol = 1, strip.position = "left")



ggsave("figures_tables/figure5_summary/annotated_proteins.png", plot = annotated_proteins, dpi = 600, width = 18 , height = 25, units = "cm", bg = "white")



####

# Summarizing amino acid data for shap and overlapping positions

overlaps <- plotting_df %>%
  group_by(segment, position) %>%
  filter(n_distinct(track) > 1 | any(track == "shap")) %>%
  ungroup() %>%
  select(segment, position, trait, track)

traits <- unique(overlaps$trait)

overlap_positions_complete <- overlaps %>%
  group_by(segment) %>%
  summarize(positions = list(unique(position))) %>%
  rowwise() %>%
  mutate(grid = list(expand.grid(
    segment = segment,
    position = positions,
    trait = traits
  ))) %>%
  unnest(cols = c(grid), names_sep = "_") %>%
  select(segment = grid_segment, position = grid_position, trait = grid_trait)


summarize_branch_diff <- function(file_path, segment_name) {
  # Read the data
  branchdiff <- read.table(file_path)
  colnames(branchdiff) <- c("node_from", "node_to", "branchlen", "position", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "trait", "external")
  
  # Calculate total residue counts
  total_residue_counts <- branchdiff %>%
    filter(branchlen <= 15) %>%
    mutate(segment = segment_name,
           position = as.numeric(position),
           trait = sub("^[^-]*-", "", trait)) %>%
    group_by(position, trait) %>%
    summarise(total_residue = n()) %>%
    ungroup()
  
  # Tally for residue_to
  residue_tally <- branchdiff %>%
    filter(branchlen <= 15) %>%
    mutate(position = as.numeric(position),
           trait = sub("^[^-]*-", "", trait)) %>%
    group_by(position, trait, residue_to) %>%
    tally() %>%
    dplyr::rename(n_residue = n) %>%
    ungroup()
  
  # Summarize residue frequencies
  residue_summary <- left_join(residue_tally, total_residue_counts, by = c("position", "trait")) %>%
    group_by(position, trait) %>%
    arrange(desc(n_residue), .by_group =  TRUE) %>%
    summarise(
      "Amino Acids (%)" = paste0(residue_to, "(", round(n_residue/total_residue*100, 1), "%)", collapse = ", ")
    ) %>%
    ungroup() %>%
    mutate(segment = segment_name, .before = "position")
  
  # Tally for mutation
  mutation_tally <- branchdiff %>%
    filter(branchlen <= 15) %>%
    mutate(position = as.numeric(position),
           trait = sub("^[^-]*-", "", trait),
           mutation = ifelse(residue_from != residue_to, paste(residue_from, "-", residue_to, sep = ""), NA)) %>%
    group_by(position, trait, mutation) %>%
    tally() %>%
    filter(!is.na(mutation)) %>%
    dplyr::rename(n_mutation = n) %>%
    ungroup()
  
  # Summarize mutation frequencies
  mutation_summary <- left_join(mutation_tally, total_residue_counts, by = c("position", "trait")) %>%
    group_by(position, trait) %>%
    arrange(desc(n_mutation), .by_group =  TRUE) %>%
    summarise(
      "Mutations (%)" = ifelse(sum(n_mutation) == 0, "NA",
                               paste0(mutation, "(", round(n_mutation/total_residue*100, 1), "%)", collapse = ", "))
    ) %>%
    ungroup() %>%
    mutate(segment = segment_name, .before = "position")
  
  # Join the summaries
  branchdiff_summarised <- left_join(residue_summary, mutation_summary, by = c("segment", "position", "trait"))
  
  return(branchdiff_summarised)
}

pb2_branchdiff_summary <- summarize_branch_diff("anclib_files/pb2/pb2_branchdiffinfo_aa_br2_all.txt", "PB2")
pb1_branchdiff_summary <- summarize_branch_diff("anclib_files/pb1/pb1_branchdiffinfo_aa_br2_all.txt", "PB1")
pb1f2_branchdiff_summary <- summarize_branch_diff("anclib_files/pb1/pb1-f2_branchdiffinfo_aa_br2_all.txt", "PB1-F2")
pa_branchdiff_summary <- summarize_branch_diff("anclib_files/pa/pa_branchdiffinfo_aa_br2_all.txt", "PA")
pax_branchdiff_summary <- summarize_branch_diff("anclib_files/pa/pa-x_branchdiffinfo_aa_br2_all.txt", "PA-X")
h1_branchdiff_summary <- summarize_branch_diff("anclib_files/ha/h1_branchdiffinfo_aa_br2_all.txt", "H1")
h3_branchdiff_summary <- summarize_branch_diff("anclib_files/ha/h3_branchdiffinfo_aa_br2_all.txt", "H3")
np_branchdiff_summary <- summarize_branch_diff("anclib_files/np/np_branchdiffinfo_aa_br2_all.txt", "NP")
n1_branchdiff_summary <- summarize_branch_diff("anclib_files/na/n1_branchdiffinfo_aa_br2_all.txt", "N1")
n2_branchdiff_summary <- summarize_branch_diff("anclib_files/na/n2_branchdiffinfo_aa_br2_all.txt", "N2")
m1_branchdiff_summary <- summarize_branch_diff("anclib_files/mp/m1_branchdiffinfo_aa_br2_all.txt", "M1")
m2_branchdiff_summary <- summarize_branch_diff("anclib_files/mp/m2_branchdiffinfo_aa_br2_all.txt", "M2")
ns1_branchdiff_summary <- summarize_branch_diff("anclib_files/ns/ns1_branchdiffinfo_aa_br2_all.txt", "NS1")
nep_branchdiff_summary <- summarize_branch_diff("anclib_files/ns/nep_branchdiffinfo_aa_br2_all.txt", "NEP")

combined_branchdiff_sumarised <- bind_rows(pb2_branchdiff_summary, pb1_branchdiff_summary, pb1f2_branchdiff_summary,
                                           pa_branchdiff_summary, pax_branchdiff_summary, h1_branchdiff_summary,
                                           h3_branchdiff_summary, np_branchdiff_summary, n1_branchdiff_summary,
                                           n2_branchdiff_summary, m1_branchdiff_summary, m2_branchdiff_summary,
                                           ns1_branchdiff_summary, nep_branchdiff_summary)

overlap_positions_amino_acids <- inner_join(combined_branchdiff_sumarised, overlap_positions_complete, by = c("segment", "position", "trait"))






