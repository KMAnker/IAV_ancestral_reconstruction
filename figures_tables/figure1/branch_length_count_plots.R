library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("PATH_TO_PROJECT_DIRECTORY")
setwd("/Users/B246311/Library/CloudStorage/OneDrive-Sundhedsdatastyrelsen/Skrivebord/IAV_ancestral_reconstruction")

trait_order <- c("swine-swine","human-human", "human-swine", "swine-human")

pb2_aa_branchdiff <- read.table("anclib_files/pb2/pb2_branchdiffinfo_aa_br2_all.txt")
colnames(pb2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pb2_summarised_branchlen <- pb2_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "PB2")

pb2_branchcounts <- pb2_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "PB2",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


pb1_aa_branchdiff <- read.table("anclib_files/pb1/pb1_branchdiffinfo_aa_br2_all.txt")
colnames(pb1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pb1_summarised_branchlen <- pb1_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "PB1")

pb1_branchcounts <- pb1_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "PB1",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


pa_aa_branchdiff <- read.table("anclib_files/pa/pa_branchdiffinfo_aa_br2_all.txt")
colnames(pa_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pa_summarised_branchlen <- pa_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "PA")

pa_branchcounts <- pa_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "PA",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


h1_aa_branchdiff <- read.table("anclib_files/ha/h1_branchdiffinfo_aa_br2_all.txt")
colnames(h1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

h1_summarised_branchlen <- h1_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "H1")

h1_branchcounts <- h1_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "H1",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


h3_aa_branchdiff <- read.table("anclib_files/ha/h3_branchdiffinfo_aa_br2_all.txt")
colnames(h3_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

h3_summarised_branchlen <- h3_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "H3")

h3_branchcounts <- h3_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "H3",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


np_aa_branchdiff <- read.table("anclib_files/np/np_branchdiffinfo_aa_br2_all.txt")
colnames(np_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

np_summarised_branchlen <- np_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "NP")

np_branchcounts <- np_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "NP",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


n1_aa_branchdiff <- read.table("anclib_files/na/n1_branchdiffinfo_aa_br2_all.txt")
colnames(n1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

n1_summarised_branchlen <- n1_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "N1")

n1_branchcounts <- n1_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "N1",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))

n2_aa_branchdiff <- read.table("anclib_files/na/n2_branchdiffinfo_aa_br2_all.txt")
colnames(n2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

n2_summarised_branchlen <- n2_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "N2")

n2_branchcounts <- n2_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "N2",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


mp_aa_branchdiff <- read.table("anclib_files/mp/m1_branchdiffinfo_aa_br2_all.txt")
colnames(mp_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

mp_summarised_branchlen <- mp_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "MP")

mp_branchcounts <- mp_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "MP",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


ns_aa_branchdiff <- read.table("anclib_files/ns/ns1_branchdiffinfo_aa_br2_all.txt")
colnames(ns_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

ns_summarised_branchlen <- ns_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-")) %>%
  select(c(trait, branchlen)) %>%
  mutate(segment = "NS")

ns_branchcounts <- ns_aa_branchdiff %>%
  mutate(trait = paste(trait_from, trait_to, sep = "-"),
         segment = "NS",
         branch_id = paste(node_from, node_to, sep = "-")) %>%
  group_by(segment, trait) %>%
  summarize(branch_count = n_distinct(branch_id),
            filtered_branch_count = n_distinct(branch_id[branchlen <= 15]))


##########################

combined_branchlength_df <- rbind(pb2_summarised_branchlen, pb1_summarised_branchlen, pa_summarised_branchlen,
                                 h1_summarised_branchlen, h3_summarised_branchlen, np_summarised_branchlen, 
                                 n1_summarised_branchlen, n2_summarised_branchlen, mp_summarised_branchlen, ns_summarised_branchlen)

combined_branchlength_df$segment <- factor(combined_branchlength_df$segment, levels= c("PB2", "PB1", "PA", "H1", "H3", "NP", "N1", "N2", "MP", "NS"))
combined_branchlength_df$trait <- factor(combined_branchlength_df$trait, levels=trait_order)

combined_branchcount_df <- rbind(pb2_branchcounts, pb1_branchcounts, pa_branchcounts,
                                 h1_branchcounts, h3_branchcounts, np_branchcounts, 
                                 n1_branchcounts, n2_branchcounts, mp_branchcounts, ns_branchcounts)

combined_branchcount_df$trait <- factor(combined_branchcount_df$trait, levels= trait_order)
combined_branchcount_df$segment <- factor(combined_branchcount_df$segment, levels= c("PB2", "PB1", "PA", "H1", "H3", "NP", "N1", "N2", "MP", "NS"))


###########################

branchlength_plot <- ggplot(combined_branchlength_df, aes(x = trait, y = branchlen, fill = trait)) +
  geom_boxplot(alpha = 1, linewidth = 0.1, color = "grey40", outlier.size = 0.1) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "grey40", size = 0.2) +
  scale_fill_manual(values = c("#45BACF", "#D78B5E", "#A5CD92", "#FFD685")) +
  labs(y = "Branch lengths") +
  facet_wrap(~segment, ncol = 5) +
  theme_classic() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 5),
        legend.spacing = unit(-0.5, "pt"),
        legend.margin = margin(unit(c(0.1,0,0.5,0), "lines")),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(linewidth = 0.5),
        strip.text = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        panel.spacing.x = unit(0.1, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))


branchcount_plot <- ggplot(combined_branchcount_df, aes(x = trait)) +
  geom_bar(aes(y = branch_count, fill = trait), stat = "identity", position = "dodge", alpha = 0.5) +
  geom_bar(aes(y = filtered_branch_count, fill = trait), stat = "identity", position = "dodge") +
  geom_text(aes(y = filtered_branch_count, label = filtered_branch_count), size = 1.5, vjust = -1) +
  labs(y = "No. of branches") +
  facet_wrap(~segment, ncol = 5) +
  scale_fill_manual(values = c("#45BACF", "#D78B5E", "#A5CD92", "#FFD685")) +
  scale_y_continuous(breaks = seq(0, 2100, by = 500), limits = c(0, 2100)) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 5),
        legend.spacing = unit(-0.5, "pt"),
        legend.margin = margin(unit(c(0.1,0,0.5,0), "lines")),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(linewidth = 0.5),
        strip.text = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        panel.spacing.x = unit(0.1, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        legend.box.spacing = unit(0.1, "lines"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  scale_x_discrete(expansion(add = c(1.5, 1.5)))


combined_branch_plots <- ggarrange(branchlength_plot, branchcount_plot, 
                                   ncol = 1, 
                                   labels = c("c", "d"), 
                                   font.label = list(size = unit(8, "pt")))


ggsave("figures_tables/figure1/branch_length_counts_plot.png", plot = combined_branch_plots, dpi = 300, width = 8 , height = 12, units = "cm" )

