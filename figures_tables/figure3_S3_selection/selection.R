library(tidyverse)
library(dplyr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)

setwd("/Users/kman/Desktop/ancestral_reconstruction_project_final/")

#######################################
### Figure 3 - absrel and meme overview

##
absrel <- read.table("results/hyphy/branch_selection.txt", header = TRUE) %>%
  mutate(Trait = case_when(
    Trait == "hu-hu" ~ "human-human",
    Trait == "sw-sw" ~ "swine-swine",
    Trait == "sw-hu" ~ "swine-human",
    Trait == "hu-sw" ~ "human-swine"))

trait_order <- c("human-human","swine-swine","human-swine", "swine-human")
protein_order <- c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")
absrel$Trait <- factor(absrel$Trait, levels=trait_order)
absrel$Segment <- factor(absrel$Segment, levels=protein_order)

absrel <- absrel %>%
  group_by(Segment) %>%
  mutate(total_selected_branches = sum(N_selected_branches),
         fraction_selected_branches_of_total = total_selected_branches/Total_branches_tree,
         normalized_fraction_selected = Fraction_selected/fraction_selected_branches_of_total)


branch_selection_overview <- ggplot(absrel, aes(x = Trait, y = normalized_fraction_selected, fill = Trait)) +
  geom_bar(stat = "identity", position = "dodge", width = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "",
       x = "",
       y = "Relative fraction of branches with positive selection") +
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  facet_wrap(~Segment, nrow=2) +
  theme_bw() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


##
meme <- read.table("results/hyphy/meme_selection.txt", header = TRUE)

trait_order <- c("human-human","swine-swine","human-swine", "swine-human")
protein_order <- c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")
meme$Trait <- factor(meme$Trait, levels=trait_order)
meme$Protein <- factor(meme$Protein, levels=protein_order)

meme_overview_plot <- meme %>%
  filter(!Protein %in% c("PB1-F2", "PA-X")) %>%
  group_by(Protein, Position) %>%
  mutate(total_branches_w_selection_at_position = sum(branches_under_selection),
         total_fraction_branches_w_selection_of_total_branches = total_branches_w_selection_at_position/Total_branches) %>%
  ungroup() %>%
  group_by(Protein,Trait) %>%
  mutate(fraction_selected_branches_of_trait_branches = branches_under_selection/Total_branches_trait,
         Relative_fraction_selected_per_trait = fraction_selected_branches_of_trait_branches/total_fraction_branches_w_selection_of_total_branches)

codon_selection_overview <- ggplot(meme_overview_plot, aes(x = Trait, y = Proportion_selected_positions, fill = Trait)) +
  geom_bar(stat = "identity", position = "dodge", width = 1) +
  labs(title = "",
       x = "",
       y = "% codons under positive selection") +
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  facet_wrap(~Protein, nrow=2) +
  theme_bw() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


###
Figure3_combined_plot <- ggarrange(branch_selection_overview, codon_selection_overview,
                               ncol = 1,
                               nrow = 2,
                               common.legend = TRUE,
                               legend = "top",
                               labels = c("a", "b"), 
                               font.label = list(size = 10))

ggsave("/Users/kman/Desktop/ancestral_reconstruction_project_final/figures_tables/figure3_S3_selection/figure3_combined_selection_overview.png", plot = Figure3_combined_plot, dpi = 300, width = 8.5 , height = 15, units = "cm" )





#######################################
### Figure S2 - meme details



# PB2

#pb2 <- meme_overview_plot %>%
#  filter(Protein == "PB2") %>%
#  group_by(Position) %>%
#  arrange(Percentage_of_trait_branches) %>%
#  ungroup()
#pb2$Trait <- factor(pb2$Trait, levels = unique(pb2$Trait))

#pb2$label <- NA
#pb2$label[pb2$Percentage_of_trait_branches > 2] <- pb2$Position[pb2$Percentage_of_trait_branches > 2]
#for (pos in unique(pb2$Position)) {
#  rows <- which(pb2$Position == pos & !is.na(pb2$label))
#  
#  if (length(rows) > 1) {
#    # Retain only the label for the bar with the highest Percentage_of_trait_branches
#    rows_to_remove <- rows[which(pb2$Percentage_of_trait_branches[rows] != max(pb2$Percentage_of_trait_branches[rows]))]
#    pb2$label[rows_to_remove] <- NA
#  }
#}

#pb2_meme_plot <- ggplot(pb2, aes(x = Position, y = Percentage_of_trait_branches, fill = Trait)) +
#  geom_bar(stat = "identity", position = "stack", width = 2) +
#  labs(title = "PB2",
#       x = "Position",
#       y = "% branches") +
#  scale_x_continuous(breaks = seq(0, 760, by = 50), limits = c(0, 760)) +
#  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
#  geom_text_repel(data = subset(pb2, Percentage_of_trait_branches > 2), 
#                  aes(label = label), position = position_stack(),  size = 2, vjust = 1, hjust = 1.1,
#                  min.segment.length = 0.1,segment.size = 0.3, max.overlaps = 20)+
#  scale_fill_manual(values = c("human-human" = "#D78B5E", 
#                               "swine-swine" = "#45BACF", 
#                               "human-swine" = "#A5CD92", 
#                               "swine-human" = "#FFD685")) +
#  theme_classic() +
#  theme(plot.title = element_text(size = 7),
#        legend.position = "top",
#        legend.title = element_blank(),
#        legend.key.size = unit(6, "points"),
#        legend.text = element_text(size = 6),
#        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
#        axis.text.y = element_text(size = 6),
#        axis.title.x = element_text(size = 7),
#        axis.title.y = element_text(size = 7),
#        axis.line = element_line(size = 0.5),
#        axis.ticks = element_line(size = 0.5),
#        strip.background = element_rect(size = 0.6),
#        panel.spacing.x = unit(0.2, "lines"), 
#        panel.spacing.y = unit(0.1, "lines"),
#        panel.grid.minor.x = element_blank(),
#        panel.grid.major.x = element_blank(),
#        legend.box.spacing = unit(0.1, "lines"),
#        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

#pb2_label_positions <- unlist(unique(na.omit(pb2$label)))
#pb2_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/pb2/pb2_mut_table.txt", header = TRUE)
#pb2_muts_at_selected_pos <- pb2_mut_data %>%
#  filter(seq_pos %in% pb2_label_positions) %>%
#  select(trait, seq_pos, mut, freq)

pb2 <- meme_overview_plot %>%
  filter(Protein == "PB2")
#pb2$Trait <- factor(pb2$Trait, levels = unique(pb2$Trait))

pb2_meme_plot <- ggplot(pb2, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "PB2",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 760, by = 50), limits = c(0, 760)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

pb2_positions <- pb2 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

pb2_mut_data <- read.table("figures_tables/distribution_plots/pb2/pb2_mut_table.txt", header = TRUE)

pb2_muts_at_selected_pos <- pb2_mut_data %>%
  inner_join(pb2_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

pb2_muts_at_selected_pos_2 <- pb2_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "PB2", .before = seq_pos)

# PB1
pb1 <- meme_overview_plot %>%
  filter(Protein == "PB1")
#pb1$Trait <- factor(pb1$Trait, levels = unique(pb1$Trait))

pb1_meme_plot <- ggplot(pb1, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "PB1",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 760, by = 50), limits = c(0, 760)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

pb1_positions <- pb1 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

pb1_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/pb1/pb1_mut_table.txt", header = TRUE)

pb1_muts_at_selected_pos <- pb1_mut_data %>%
  inner_join(pb1_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

pb1_muts_at_selected_pos_2 <- pb1_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "PB1", .before = seq_pos)


# PA
pa <- meme_overview_plot %>%
  filter(Protein == "PA")
#pa$Trait <- factor(pa$Trait, levels = unique(pa$Trait))

pa_meme_plot <- ggplot(pa, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "PA",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 760, by = 50), limits = c(0, 760)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

pa_positions <- pa %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

pa_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/pa/pa_mut_table.txt", header = TRUE)

pa_muts_at_selected_pos <- pa_mut_data %>%
  inner_join(pa_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

pa_muts_at_selected_pos_2 <- pa_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "PA", .before = seq_pos)

# H1
h1 <- meme_overview_plot %>%
  filter(Protein == "H1")
#h1$Trait <- factor(h1$Trait, levels = unique(h1$Trait))

h1_meme_plot <- ggplot(h1, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "H1",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 570, by = 50), limits = c(0, 570)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, margin = margin(unit(c(0,12,0,0), "lines"))),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

h1_positions <- h1 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

h1_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/ha/h1_mut_table.txt", header = TRUE)

h1_muts_at_selected_pos <- h1_mut_data %>%
  inner_join(h1_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

h1_muts_at_selected_pos_2 <- h1_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "H1", .before = seq_pos)


# H3
h3 <- meme_overview_plot %>%
  filter(Protein == "H3")
#h3$Trait <- factor(h3$Trait, levels = unique(h3$Trait))

h3_meme_plot <- ggplot(h3, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "H3",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 570, by = 50), limits = c(0, 570)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

h3_positions <- h3 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

h3_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/ha/h3_mut_table.txt", header = TRUE)

h3_muts_at_selected_pos <- h3_mut_data %>%
  inner_join(h3_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

h3_muts_at_selected_pos_2 <- h3_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "H3", .before = seq_pos)

# NP
np <- meme_overview_plot %>%
  filter(Protein == "NP")
#np$Trait <- factor(np$Trait, levels = unique(np$Trait))

np_meme_plot <- ggplot(np, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "NP",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 510, by = 50), limits = c(0, 510)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

np_positions <- np %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

np_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/np/np_mut_table.txt", header = TRUE)

np_muts_at_selected_pos <- np_mut_data %>%
  inner_join(np_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

np_muts_at_selected_pos_2 <- np_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "NP", .before = seq_pos)

# N1
n1 <- meme_overview_plot %>%
  filter(Protein == "N1")
#n1$Trait <- factor(n1$Trait, levels = unique(n1$Trait))

n1_meme_plot <- ggplot(n1, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "N1",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 470, by = 50), limits = c(0, 470)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

n1_positions <- n1 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

n1_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/na/n1_mut_table.txt", header = TRUE)

n1_muts_at_selected_pos <- n1_mut_data %>%
  inner_join(n1_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

n1_muts_at_selected_pos_2 <- n1_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "N1", .before = seq_pos)

# N2
n2 <- meme_overview_plot %>%
  filter(Protein == "N2")
#n2$Trait <- factor(n2$Trait, levels = unique(n2$Trait))

n2_meme_plot <- ggplot(n2, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "N2",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 470, by = 50), limits = c(0, 470)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

n2_positions <- n2 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

n2_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/na/n2_mut_table.txt", header = TRUE)

n2_muts_at_selected_pos <- n2_mut_data %>%
  inner_join(n2_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

n2_muts_at_selected_pos_2 <- n2_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "N2", .before = seq_pos)

# M1
m1 <- meme_overview_plot %>%
  filter(Protein == "M1")
#m1$Trait <- factor(m1$Trait, levels = unique(m1$Trait))

m1_meme_plot <- ggplot(m1, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "M1",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 260, by = 50), limits = c(0, 260)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

m1_positions <- m1 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

m1_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/mp/m1_mut_table.txt", header = TRUE)

m1_muts_at_selected_pos <- m1_mut_data %>%
  inner_join(m1_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

m1_muts_at_selected_pos_2 <- m1_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "M1", .before = seq_pos)

# M2
m2 <- meme_overview_plot %>%
  filter(Protein == "M2")
#m2$Trait <- factor(m2$Trait, levels = unique(m2$Trait))

m2_meme_plot <- ggplot(m2, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "M2",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 100, by = 50), limits = c(0, 100)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

m2_positions <- m2 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

m2_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/mp/m2_mut_table.txt", header = TRUE)

m2_muts_at_selected_pos <- m2_mut_data %>%
  inner_join(m2_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

m2_muts_at_selected_pos_2 <- m2_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "M2", .before = seq_pos)

# NS1
ns1 <- meme_overview_plot %>%
  filter(Protein == "NS1")
#ns1$Trait <- factor(ns1$Trait, levels = unique(ns1$Trait))

ns1_meme_plot <- ggplot(ns1, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "NS1",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 240, by = 50), limits = c(0, 240)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

ns1_positions <- ns1 %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

ns1_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/ns/ns1_mut_table.txt", header = TRUE)

ns1_muts_at_selected_pos <- ns1_mut_data %>%
  inner_join(ns1_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

ns1_muts_at_selected_pos_2 <- ns1_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "NS1", .before = seq_pos)

# NEP
nep <- meme_overview_plot %>%
  filter(Protein == "NEP")
#nep$Trait <- factor(nep$Trait, levels = unique(nep$Trait))

nep_meme_plot <- ggplot(nep, aes(x = Position, y = Relative_fraction_selected_per_trait, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", width = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "NEP",
       x = "Position",
       y = "Relative fraction of branches with selection at site for each trait") +
  scale_x_continuous(breaks = seq(0, 130, by = 50), limits = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  scale_fill_manual(values = c("human-human" = "#D78B5E", 
                               "swine-swine" = "#45BACF", 
                               "human-swine" = "#A5CD92", 
                               "swine-human" = "#FFD685")) +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle=45, hjust = 0.9),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))

nep_positions <- nep %>%
  filter(Relative_fraction_selected_per_trait > 1) %>%
  select(Position, Trait, Relative_fraction_selected_per_trait)

nep_mut_data <- read.table("figures_tables/distribution_blosum_plots/distribution_plots/ns/nep_mut_table.txt", header = TRUE)

nep_muts_at_selected_pos <- nep_mut_data %>%
  inner_join(nep_positions, by = c("seq_pos" = "Position", "trait" = "Trait")) %>%
  select(trait, seq_pos, mut, freq, Relative_fraction_selected_per_trait)

nep_muts_at_selected_pos_2 <- nep_muts_at_selected_pos %>%
  group_by(seq_pos, trait) %>%
  arrange(desc(freq), .by_group =  TRUE) %>%
  summarise(Mutations = paste(mut, "(", freq, ")", sep = "", collapse = ", ")) %>%
  ungroup()%>%
  mutate(Segment = "NEP", .before = seq_pos)

################################################################################

combined_meme_dist_1 <- ggarrange(pb2_meme_plot, pb1_meme_plot, pa_meme_plot, 
                               ncol = 2,
                               nrow = 2,
                               common.legend = TRUE,
                               legend = "top") +
  theme(plot.margin = unit(c(0,0,0,0.1), "null")) 


combined_meme_dist_2 <- ggarrange(h1_meme_plot, h3_meme_plot,
                             ncol = 2,
                             nrow = 1,
                             widths = c(1.1, 0.9)) +
  theme(plot.margin = unit(c(0,0.3,0,0.058), "null")) 


combined_meme_dist_3 <- ggarrange(np_meme_plot, n1_meme_plot, n2_meme_plot,
                             ncol = 3,
                             nrow = 1) +
  theme(plot.margin = unit(c(0,0,0,0.1), "null"))  

combined_meme_dist_4 <- ggarrange(m1_meme_plot, m2_meme_plot, ns1_meme_plot, nep_meme_plot, 
                             ncol = 4,
                             nrow = 1,
                             widths = c(0.32,0.2,0.31,0.22)) +
  theme(plot.margin = unit(c(0,0.4,0,0.14), "null")) 


combined_meme_dist_all <- ggarrange(combined_meme_dist_1, combined_meme_dist_2, combined_meme_dist_3, combined_meme_dist_4,
                               ncol = 1,
                               nrow = 4,
                               common.legend = TRUE,
                               legend = "top",
                               heights = c(2,1,1,1))


ggsave("/Users/kman/Desktop/ancestral_reconstruction_project_final/figures_tables/figure3_S3_selection/selected_positions_distribution.png", plot = combined_meme_dist_all, dpi = 300, width = 16 , height = 20, units = "cm" )



combined_mut_list <- rbind(pb2_muts_at_selected_pos_2, pb1_muts_at_selected_pos_2, pa_muts_at_selected_pos_2,
                           h1_muts_at_selected_pos_2, h3_muts_at_selected_pos_2, np_muts_at_selected_pos_2,
                           n1_muts_at_selected_pos_2, n2_muts_at_selected_pos_2, m1_muts_at_selected_pos_2,
                           m2_muts_at_selected_pos_2, ns1_muts_at_selected_pos_2, nep_muts_at_selected_pos_2)

write.table(combined_mut_list, "/Users/kman/Desktop/ancestral_reconstruction_project_final/figures_tables/figure3_S3_selection/selected_positions_mut_list.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


################################################################################

# PB1-F2
pb1f2 <- meme %>%
  filter(Protein == "PB1-F2")


pb1f2_meme_plot <- ggplot(pb1f2, aes(x = Position, y = Percentage_of_trait_branches, fill = Trait)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.2, width = 0.5), width = 1.5) +
  labs(title = "PB1-F2",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(0, 101, by = 10), limits = c(0, 101)) +
  geom_text(data = subset(pb1f2, Percentage_of_trait_branches > 1), 
            aes(label = Position, y = Percentage_of_trait_branches), 
            position = position_dodge2(preserve = "total", width = 0.5), size = 2, vjust = -0.1, hjust = 0.1, angle = 45)+
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        axis.text.x = element_text(size = 6, vjust = 0.8, angle = 45),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


# PA-X
pax <- meme %>%
  filter(Protein == "PA-X")

pax_meme_plot <- ggplot(pax, aes(x = Position, y = Percentage_of_trait_branches, fill = Trait)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.2, width = 0.5), width = 3) +
  labs(title = "PA-X",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(0, 260, by = 50), limits = c(1, 260)) +
  geom_text(data = subset(pax, Percentage_of_trait_branches > 1), 
            aes(label = Position, y = Percentage_of_trait_branches), 
            position = position_dodge2(preserve = "total", width = 0.1), size = 2, vjust = -0.5, hjust = -0.5, angle = 45)+
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  theme_classic() +
  theme(legend.position = "",
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        axis.text.x = element_text(size = 7, vjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))



########

meme2 <- read.table("results/hyphy/meme_selection_2.txt", header = TRUE)

trait_order <- c("hu-hu","sw-sw","hu-sw", "sw-hu")
protein_order <- c("PB2", "PB1", "PA", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")
meme2$trait <- factor(meme2$trait, levels=trait_order)
meme2$Segment <- factor(meme2$Segment, levels=protein_order)


# Count how many traits each position is associated with within each segment
position_counts <- meme2 %>%
  group_by(Segment, Position) %>%
  summarise(NumTraits = n_distinct(trait), .groups = 'drop')
# Merge back with the original data and categorize positions as unique or shared
categorized_data <- meme2 %>%
  left_join(position_counts, by = c("Segment", "Position")) %>%
  mutate(Status = ifelse(NumTraits > 1, "Shared", "Unique"))
# Summarize the count of unique and shared positions for each trait within each segment
trait_summary <- categorized_data %>%
  group_by(Segment, trait, Status) %>%
  summarise(Count = n_distinct(Position), .groups = 'drop')

proportion_shared <- ggplot(trait_summary, aes(x = trait, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Segment, nrow = 2) +
  labs(y = "Proportion of selected positions", x = "Trait", fill = "Position Type") +
  scale_fill_manual(values = c("#B7B7A6", "#91A178")) +
  theme_bw() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1.1),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


# Identify shared positions and associated traits
shared_positions <- meme2 %>%
  group_by(Segment, Position) %>%
  summarise(Traits = toString(sort(unique(trait))), .groups = 'drop') %>%
  filter(n_distinct(Traits) > 1)
# Count occurrences of each shared trait combination
shared_counts <- shared_positions %>%
  group_by(Segment, Traits) %>%
  summarise(Count = n(), .groups = 'drop')

shared_counts$Traits <- factor(shared_counts$Traits, levels=c("hu-hu", "sw-sw", "hu-sw", "sw-hu", 
                                                              "hu-hu, sw-sw", "hu-hu, hu-sw", "hu-hu, sw-hu",
                                                              "sw-sw, hu-sw", "sw-sw, sw-hu",
                                                              "hu-sw, sw-hu",
                                                              "sw-sw, hu-sw, sw-hu"))
distribution_shared <- ggplot(shared_counts, aes(x = Traits, y = Count, fill = Traits)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Segment, nrow = 2) +
  labs(y = "Count of selected positions", fill = "Transmission types") +
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685", "#8EA296", "#BEAC78", "#EBB071", "#75C3B0", 
                               "#A2C8AA", "#D2D18B", "#A3C3C1")) +
  theme_bw() +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(6, "points"),
        legend.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


combined_shared_positions <- ggarrange(proportion_shared, distribution_shared,
                               ncol = 1,
                               nrow = 2,
                               labels = c("a", "b"), 
                               font.label = list(size = 10))

combined_selection_positions_figure <- ggarrange(combined_shared_positions, combined_meme_dist_all,
                                       ncol = 2,
                                       nrow = 1,
                                       labels = c("a", "c"),
                                       widths = c(0.85, 1.15),
                                       font.label = list(size = 10))

ggsave("/Users/kman/Desktop/ancestral_reconstruction_project_final/figures_tables/figure3_S3_selection/figS3_combined_selection_positions_overview.png", plot = combined_selection_positions_figure, dpi = 300, width = 17 , height = 20, units = "cm" )

