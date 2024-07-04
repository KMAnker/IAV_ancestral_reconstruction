library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("PATH_TO_PROJECT_DIRECTORY")

#######################################
### Figure 3 - absrel and meme overview

##
absrel <- read.table("selection_analysis/branch_selection.txt", header = TRUE) %>%
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
  theme_classic() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        legend.text = element_text(size = 7),
        legend.spacing = unit(-0.5, "pt"),
        legend.margin = margin(unit(c(0.1,0,0.5,0.1), "lines")),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.2),
        strip.background = element_rect(linewidth = 0.5),
        strip.text = element_text(size = 8),
        panel.spacing.x = unit(0.1, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))


##
meme <- read.table("selection_analysis/meme_selection.txt", header = TRUE)

trait_order <- c("human-human","swine-swine","human-swine", "swine-human")
protein_order <- c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")
meme$trait <- factor(meme$trait, levels=trait_order)
meme$protein <- factor(meme$protein, levels=protein_order)

meme <- meme %>%
  filter(!protein == "PB1-F2")

codon_selection_overview <- ggplot(meme, aes(x = trait, y = proportion_selected_positions, fill = trait)) +
  geom_bar(stat = "identity", position = "dodge", width = 1) +
  labs(title = "",
       x = "",
       y = "% codons under positive selection") +
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685")) +
  facet_wrap(~protein, nrow=2) +
  theme_classic() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        legend.text = element_text(size = 7),
        legend.spacing = unit(-0.5, "pt"),
        legend.margin = margin(unit(c(0.1,0,0.5,0.1), "lines")),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.2),
        strip.background = element_rect(linewidth = 0.5),
        strip.text = element_text(size = 8),
        panel.spacing.x = unit(0.1, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))

###
Figure3_combined_plot <- ggarrange(branch_selection_overview, codon_selection_overview,
                               ncol = 1,
                               nrow = 2,
                               common.legend = TRUE,
                               legend = "top",
                               labels = c("a", "b"), 
                               font.label = list(size = 10))

ggsave("figures_tables/figure3_S3S4_selection/figure3_combined_selection_overview.png", plot = Figure3_combined_plot, dpi = 300, width = 8.5 , height = 15, units = "cm" )


################################################################################

# Sup. figures - unique vs shared selected positions

meme2 <- meme %>%
  select(protein, position, trait) %>%
  mutate(trait = case_when(
    trait == "human-human" ~ "hu-hu",
    trait == "swine-swine" ~ "sw-sw",
    trait == "swine-human" ~ "sw-hu",
    trait == "human-swine" ~ "hu-sw"))

trait_order <- c("hu-hu","sw-sw","hu-sw", "sw-hu")
protein_order <- c("PB2", "PB1", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")
meme2$trait <- factor(meme2$trait, levels=trait_order)
meme2$protein <- factor(meme2$protein, levels=protein_order)

# Count how many traits each position is associated with within each segment
position_counts <- meme2 %>%
  group_by(protein, position) %>%
  summarise(num_traits = n_distinct(trait), .groups = 'drop')
# Merge back with the original data and categorize positions as unique or shared
categorized_data <- meme2 %>%
  left_join(position_counts, by = c("protein", "position")) %>%
  mutate(status = ifelse(num_traits > 1, "shared", "unique"))
# Summarize the count of unique and shared positions for each trait within each segment
trait_summary <- categorized_data %>%
  group_by(protein, trait, status) %>%
  summarise(count = n_distinct(position), .groups = 'drop')

proportion_shared <- ggplot(trait_summary, aes(x = trait, y = count, fill = status)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~protein, ncol = 4) +
  labs(y = "Proportion of selected positions", x = "Trait") +
  scale_fill_manual(values = c("#B7B7A6", "#91A178")) +
  theme_classic() +
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
        panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


# Identify shared positions and associated traits
shared_positions <- meme2 %>%
  group_by(protein, position) %>%
  summarise(traits = toString(sort(unique(trait))), .groups = 'drop') %>%
  filter(n_distinct(traits) > 1) %>%
  mutate(traits = str_replace_all(traits, ",", ";"))
# Count occurrences of each shared trait combination
shared_counts <- shared_positions %>%
  group_by(protein, traits) %>%
  summarise(count = n(), .groups = 'drop')

shared_counts$traits <- factor(shared_counts$traits, levels=c("hu-hu", "sw-sw", "hu-sw", "sw-hu", 
                                                              "hu-hu; sw-sw", "hu-hu; hu-sw", "hu-hu; sw-hu",
                                                              "sw-sw; hu-sw", "sw-sw; sw-hu",
                                                              "hu-hu; sw-sw; hu-sw",
                                                              "hu-hu; sw-sw; sw-hu",
                                                              "sw-sw; hu-sw; sw-hu",
                                                              "hu-hu; sw-sw; hu-sw; sw-hu"))

distribution_shared <- ggplot(shared_counts, aes(x = traits, y = count, fill = traits)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~protein, ncol = 2) +
  labs(y = "Count of selected positions", fill = "Transmission types") +
  scale_fill_manual(values = c("#D78B5E", "#45BACF", "#A5CD92", "#FFD685", 
                               "#8EA296", "#B1A175", "#EBB071", 
                               "#75C3B0", "#A2C8AA", 
                               "#91A885", 
                               "#D2D18B", 
                               "#A3C3C1", 
                               "#B0AF70")) +
  theme_classic() +
  guides(fill=guide_legend(nrow=4))+
  theme(plot.title = element_blank(),
        legend.position = "top",
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
        panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.2,0,0.2), "lines"))


combined_shared_proportions <- ggarrange(proportion_shared, distribution_shared,
                               ncol = 2,
                               nrow = 1,
                               labels = c("a", "b"),
                               widths = c(1,1.5),
                               font.label = list(size = 10))



ggsave("figures_tables/figure3_S3S4_selection/figureS4_selection_shared_overview.png", plot = combined_shared_proportions, dpi = 300, width = 17 , height = 20, units = "cm" )


#######################################

# Table of positively selected positions and mutations 

mut_data <- read.table("figures_tables/figure2_mut_count_plots/mut_table.txt", header = TRUE)

muts_at_selected_pos <- meme %>%
  select(protein, position, trait) %>%
  inner_join(mut_data, by = c("protein" = "protein", "position" = "seq_pos", "trait" = "trait")) %>%
  select(protein, trait, position, mut, mutfreq)

selection_mut_list <- muts_at_selected_pos %>%
  group_by(protein, position, trait) %>%
  arrange(desc(mutfreq), .by_group =  TRUE) %>%
  summarise(mutations = paste(mut, "(", mutfreq, ")", sep = "", collapse = ", ")) %>%
  ungroup() %>%
  mutate(trait = case_when(
    trait == "human-human" ~ "hu-hu",
    trait == "swine-swine" ~ "sw-sw",
    trait == "swine-human" ~ "sw-hu",
    trait == "human-swine" ~ "hu-sw"))

protein_order <- c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")

selection_mut_list$protein <- factor(selection_mut_list$protein, levels=protein_order)
selection_mut_list <- selection_mut_list %>%
  arrange(protein)

write.table(selection_mut_list, "figures_tables/figure3_S3_selection/selected_positions_mut_list.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


