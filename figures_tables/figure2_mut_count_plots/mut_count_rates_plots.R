library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Biostrings)
library(plotrix)

setwd("PATH_TO_PROJECT_DIRECTORY")

trait_order <- c("swine-swine","human-human", "human-swine", "swine-human")
protein_order <- c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")

theme1 <- theme(legend.position = "top",
                legend.title = element_blank(),
                legend.key.size = unit(7, "pt"),
                legend.text = element_text(size = 7),
                legend.spacing = unit(-0.5, "pt"),
                legend.margin = margin(unit(c(0,0,10,0), "lines")),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=7),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 7),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                strip.background = element_rect(linewidth = 0.5),
                strip.text = element_text(size = 8),
                axis.line = element_line(size = 0.1),
                panel.spacing.x = unit(0.1, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                legend.box.spacing = unit(0.1, "lines"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))

theme2 <- theme(axis.text.x = element_text(size = 7, vjust = 0.8, angle = 45, hjust = 0.9),
                axis.text.y = element_text(size=7),
                axis.title.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.ticks.y = element_blank(),
                strip.background = element_rect(linewidth = 0.5),
                strip.text = element_text(size = 7),
                axis.line = element_line(size = 0.1),
                panel.spacing.x = unit(0.1, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                legend.box.spacing = unit(0.1, "lines"),
                plot.margin = unit(c(0.5,2.5,0.5,2.5), "lines"))

theme3 <- theme(legend.position = "top",
                legend.title = element_blank(),
                legend.key.size = unit(7, "pt"),
                legend.text = element_text(size = 7),
                legend.spacing = unit(-0.5, "pt"),
                legend.margin = margin(unit(c(0.5,0,5,0), "lines")),
                axis.text.x = element_text(size = 7, vjust = 0.8, angle = 45, hjust = 0.9),
                axis.text.y = element_text(size=7),
                axis.title.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.ticks.y = element_blank(),
                strip.background = element_rect(linewidth = 0.5),
                strip.text = element_text(size = 7),
                axis.line = element_line(size = 0.1),
                panel.spacing.x = unit(0.1, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                legend.box.spacing = unit(0.1, "lines"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))

################################################################################

## loading data

#pb2
pb2_aa_branchdiff <- read.table("anclib_files/pb2/pb2_branchdiffinfo_aa_br2.txt")
colnames(pb2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
pb2_aa_branchdiff_filtered <- pb2_aa_branchdiff %>%
  filter(branchlen <= 15)

pb2_total_summary <- pb2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

pb2_summarised <- pb2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "PB2",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = pb2_total_summary$average_mut_total,
    std_error_total = pb2_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#pb1
pb1_aa_branchdiff <- read.table("anclib_files/pb1/pb1_branchdiffinfo_aa_br2.txt")
colnames(pb1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
pb1_aa_branchdiff_filtered <- pb1_aa_branchdiff %>%
  filter(branchlen <= 15)

pb1_total_summary <- pb1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

pb1_summarised <- pb1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "PB1",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = pb1_total_summary$average_mut_total,
    std_error_total = pb1_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#pb1-f2
pb1f2_aa_branchdiff <- read.table("anclib_files/pb1/pb1-f2_branchdiffinfo_aa_br2.txt")
colnames(pb1f2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
pb1f2_aa_branchdiff_filtered <- pb1f2_aa_branchdiff %>%
  filter(branchlen <= 15)

pb1f2_total_summary <- pb1f2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

pb1f2_summarised <- pb1f2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "PB1-F2",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = pb1f2_total_summary$average_mut_total,
    std_error_total = pb1f2_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#pa
pa_aa_branchdiff <- read.table("anclib_files/pa/pa_branchdiffinfo_aa_br2.txt")
colnames(pa_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
pa_aa_branchdiff_filtered <- pa_aa_branchdiff %>%
  filter(branchlen <= 15)

pa_total_summary <- pa_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

pa_summarised <- pa_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "PA",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = pa_total_summary$average_mut_total,
    std_error_total = pa_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#pa-x
pax_aa_branchdiff <- read.table("anclib_files/pa/pa-x_branchdiffinfo_aa_br2.txt")
colnames(pax_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
pax_aa_branchdiff_filtered <- pax_aa_branchdiff %>%
  filter(branchlen <= 15)

pax_total_summary <- pax_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

pax_summarised <- pax_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "PA-X",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = pax_total_summary$average_mut_total,
    std_error_total = pax_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#h1
h1_aa_branchdiff <- read.table("anclib_files/ha/h1_branchdiffinfo_aa_br2.txt")
colnames(h1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
h1_aa_branchdiff_filtered <- h1_aa_branchdiff %>%
  filter(branchlen <= 15)

h1_total_summary <- h1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

h1_summarised <- h1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "H1",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = h1_total_summary$average_mut_total,
    std_error_total = h1_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#h3
h3_aa_branchdiff <- read.table("anclib_files/ha/h3_branchdiffinfo_aa_br2.txt")
colnames(h3_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
h3_aa_branchdiff_filtered <- h3_aa_branchdiff %>%
  filter(branchlen <= 15)

h3_total_summary <- h3_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

h3_summarised <- h3_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "H3",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = h3_total_summary$average_mut_total,
    std_error_total = h3_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#np
np_aa_branchdiff <- read.table("anclib_files/np/np_branchdiffinfo_aa_br2.txt")
colnames(np_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
np_aa_branchdiff_filtered <- np_aa_branchdiff %>%
  filter(branchlen <= 15)

np_total_summary <- np_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

np_summarised <- np_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "NP",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = np_total_summary$average_mut_total,
    std_error_total = np_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#n1
n1_aa_branchdiff <- read.table("anclib_files/na/n1_branchdiffinfo_aa_br2.txt")
colnames(n1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
n1_aa_branchdiff_filtered <- n1_aa_branchdiff %>%
  filter(branchlen <= 15)

n1_total_summary <- n1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

n1_summarised <- n1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "N1",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = n1_total_summary$average_mut_total,
    std_error_total = n1_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#n2
n2_aa_branchdiff <- read.table("anclib_files/na/n2_branchdiffinfo_aa_br2.txt")
colnames(n2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
n2_aa_branchdiff_filtered <- n2_aa_branchdiff %>%
  filter(branchlen <= 15)

n2_total_summary <- n2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

n2_summarised <- n2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "N2",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = n2_total_summary$average_mut_total,
    std_error_total = n2_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#m1
m1_aa_branchdiff <- read.table("anclib_files/mp/m1_branchdiffinfo_aa_br2.txt")
colnames(m1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
m1_aa_branchdiff_filtered <- m1_aa_branchdiff %>%
  filter(branchlen <= 15)

m1_total_summary <- m1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

m1_summarised <- m1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "M1",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = m1_total_summary$average_mut_total,
    std_error_total = m1_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#m2
m2_aa_branchdiff <- read.table("anclib_files/mp/m2_branchdiffinfo_aa_br2.txt")
colnames(m2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
m2_aa_branchdiff_filtered <- m2_aa_branchdiff %>%
  filter(branchlen <= 15)

m2_total_summary <- m2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

m2_summarised <- m2_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "M2",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = m2_total_summary$average_mut_total,
    std_error_total = m2_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#ns1
ns1_aa_branchdiff <- read.table("anclib_files/ns/ns1_branchdiffinfo_aa_br2.txt")
colnames(ns1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
ns1_aa_branchdiff_filtered <- ns1_aa_branchdiff %>%
  filter(branchlen <= 15)

ns1_total_summary <- ns1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

ns1_summarised <- ns1_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "NS1",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = ns1_total_summary$average_mut_total,
    std_error_total = ns1_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)


#nep
nep_aa_branchdiff <- read.table("anclib_files/ns/nep_branchdiffinfo_aa_br2.txt")
colnames(nep_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")
nep_aa_branchdiff_filtered <- nep_aa_branchdiff %>%
  filter(branchlen <= 15)

nep_total_summary <- nep_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  group_by(branch_id = paste(node_from, node_to, sep = "-")) %>%
  summarize(
    mutations_per_branch = sum(residue_from != residue_to)
  ) %>%
  ungroup() %>%
  summarize(
    total_mutations = sum(mutations_per_branch),
    total_branch_count = n(),
    average_mut_total = mean(mutations_per_branch),
    std_error_total = std.error(mutations_per_branch)
  )

nep_summarised <- nep_aa_branchdiff_filtered %>%
  filter(residue_from != 'X', residue_to != 'X') %>%
  mutate(
    trait = paste(trait_from, trait_to, sep = "-"),
    protein = "NEP",
    branch_id = paste(node_from, node_to, sep = "-")
  ) %>%
  group_by(trait, protein, branch_id) %>%
  summarize(
    num_mutations = sum(residue_from != residue_to)) %>%
  group_by(trait, protein) %>%
  summarize(
    average_mut_trait = mean(num_mutations),
    std_error = std.error(num_mutations)) %>%
  mutate(
    average_mut_total = nep_total_summary$average_mut_total,
    std_error_total = nep_total_summary$std_error_total,
    normalized_mut_per_trait = average_mut_trait/average_mut_total,
    normalized_se = std_error/average_mut_total)



#####

combined_summarised <- rbind(pb2_summarised, pb1_summarised, pb1f2_summarised, pa_summarised, pax_summarised,
                             h1_summarised, h3_summarised, np_summarised, n1_summarised, n2_summarised,
                             m1_summarised, m2_summarised, ns1_summarised, nep_summarised)

combined_summarised$trait <- factor(combined_summarised$trait, levels=trait_order)
combined_summarised$protein <- factor(combined_summarised$protein, levels=protein_order)

normalized <- ggplot(combined_summarised, aes(x = trait, y = normalized_mut_per_trait, fill = trait)) +
  geom_col() +
  geom_errorbar(aes(ymin = normalized_mut_per_trait - normalized_se, ymax = normalized_mut_per_trait + normalized_se),
                width = 0.3, linewidth = 0.5, color = "grey40") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 0.3) +
  scale_fill_manual(name = "Transmission type", values = c("#45BACF", "#D78B5E", "#A5CD92", "#FFD685")) +
  labs(y = "Normalized Mutation Count") +
  facet_wrap(~protein, ncol = 7) +
  theme_classic()+ 
  theme1
normalized <- normalized + guides(fill = guide_legend(nrow = 2, byrow = TRUE))


combined_summarised_2 <- combined_summarised %>%
  group_by(protein) %>%
  summarize(
    average_mut_total = mean(average_mut_total),
    std_error_total = mean(std_error_total))

average <- ggplot(combined_summarised_2, aes(x = protein, y = average_mut_total)) +
  geom_col(fill = "#91A178") +
  geom_errorbar(aes(ymin = average_mut_total - std_error_total, ymax = average_mut_total + std_error_total), width = 0.2, size = 0.5, color = "grey40") +
  labs(y = "Average Mutation Count", x = "") +
  theme_classic()+ 
  theme2

total_and_normalized_mut_count_plots <- ggarrange(average, normalized, 
                           ncol = 1,
                           heights = c(0.6,1.4),
                           labels = c("a", "b"), 
                           font.label = list(size = 10))


################################################################################
## Bayesian analysis ##

# Read in tibbles

pb2_tibble = readRDS("bayesian_mut_analysis/output_files/pb2_comparison.rds") %>%
  mutate(segment = "PB2")
pb1_tibble = readRDS("bayesian_mut_analysis/output_files/pb1_comparison.rds") %>%
  mutate(segment = "PB1")
pb1_f2_tibble = readRDS("bayesian_mut_analysis/output_files/pb1-f2_comparison.rds") %>%
  mutate(segment = "PB1-F2")
pa_tibble = readRDS("bayesian_mut_analysis/output_files/pa_comparison.rds") %>%
  mutate(segment = "PA")
pa_x_tibble = readRDS("bayesian_mut_analysis/output_files/pa-x_comparison.rds") %>%
  mutate(segment = "PA-X")
h1_tibble = readRDS("bayesian_mut_analysis/output_files/h1_comparison.rds") %>%
  mutate(segment = "H1")
h3_tibble = readRDS("bayesian_mut_analysis/output_files/h3_comparison.rds") %>%
  mutate(segment = "H3")
np_tibble = readRDS("bayesian_mut_analysis/output_files/np_comparison.rds") %>%
  mutate(segment = "NP")
n1_tibble = readRDS("bayesian_mut_analysis/output_files/n1_comparison.rds") %>%
  mutate(segment = "N1")
n2_tibble = readRDS("bayesian_mut_analysis/output_files/n2_comparison.rds") %>%
  mutate(segment = "N2")
m1_tibble = readRDS("bayesian_mut_analysis/output_files/m1_comparison.rds") %>%
  mutate(segment = "M1")
m2_tibble = readRDS("bayesian_mut_analysis/output_files/m2_comparison.rds") %>%
  mutate(segment = "M2")
ns1_tibble = readRDS("bayesian_mut_analysis/output_files/ns1_comparison.rds") %>%
  mutate(segment = "NS1")
nep_tibble = readRDS("bayesian_mut_analysis/output_files/nep_comparison.rds") %>%
  mutate(segment = "NEP")

combined_tibble <- bind_rows(pb2_tibble, pb1_tibble, pb1_f2_tibble, pa_tibble, pa_x_tibble,
                             h1_tibble, h3_tibble, np_tibble, n1_tibble, n2_tibble, 
                             m1_tibble, m2_tibble, ns1_tibble, nep_tibble)


combined_data <- bind_rows(
  combined_tibble %>% 
    filter(class1 == "swsw", class2 == "husw", P_2_larger_1 > 0.95) %>% 
    mutate(comparison = "husw_vs_swsw"),
  combined_tibble %>% 
    filter(class1 == "huhu", class2 == "swhu", P_2_larger_1 > 0.95) %>%
    mutate(comparison = "swhu_vs_huhu"),
  combined_tibble %>% 
    filter(class1 == "swsw", class2 == "swhu", P_2_larger_1 > 0.95) %>% 
    mutate(comparison = "swhu_vs_swsw"),
  combined_tibble %>% 
    filter(class1 == "huhu", class2 == "husw", P_2_larger_1 > 0.95) %>% 
    mutate(comparison = "husw_vs_huhu"),
  combined_tibble %>% 
    filter(class1 == "huhu", class2 == "swsw", P_1_larger_2 > 0.95) %>% 
    mutate(comparison = "huhu_vs_swsw"),
  combined_tibble %>% 
    filter(class1 == "huhu", class2 == "swsw", P_2_larger_1 > 0.95) %>% 
    mutate(comparison = "swsw_vs_huhu"),
  combined_tibble %>% 
    filter(class1 == "husw", class2 == "swhu", P_2_larger_1 > 0.95) %>% 
    mutate(comparison = "swhu_vs_husw"),
  combined_tibble %>% 
    filter(class1 == "husw", class2 == "swhu", P_1_larger_2 > 0.95) %>% 
    mutate(comparison = "husw_vs_swhu"),
  combined_tibble %>% 
    filter(class1 == "huhu", class2 == "husw", P_1_larger_2 > 0.95) %>% 
    mutate(comparison = "huhu_vs_husw"),
  combined_tibble %>% 
    filter(class1 == "huhu", class2 == "swhu", P_1_larger_2 > 0.95) %>% 
    mutate(comparison = "huhu_vs_swhu"),
  combined_tibble %>% 
    filter(class1 == "swsw", class2 == "husw", P_1_larger_2 > 0.95) %>% 
    mutate(comparison = "swsw_vs_husw"),
  combined_tibble %>%
    filter(class1 == "swsw", class2 == "swhu", P_1_larger_2 > 0.95) %>% 
    mutate(comparison = "swsw_vs_swhu")
)

protein_order <- c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP")
combined_data$segment <- factor(combined_data$segment, levels=protein_order)

summary <- combined_data %>%
  group_by(comparison, segment) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(comparison, segment, fill = list(count = 0))


cross_host_comparisons <- c(
  "husw_vs_swsw", "swsw_vs_husw",
  "swhu_vs_huhu", "huhu_vs_swhu",
  "husw_vs_huhu", "huhu_vs_husw",
  "swhu_vs_swsw", "swsw_vs_swhu")

within_host_comparisons <- c("huhu_vs_swsw", "swsw_vs_huhu")


cross_host_data <- summary %>%
  filter(comparison %in% cross_host_comparisons) %>%
  mutate(
    # Relabel the comparison values
    comparison = case_when(
      comparison == "husw_vs_swsw" ~ "human-swine > swine-swine",
      comparison == "swsw_vs_husw" ~ "human-swine < swine-swine",
      comparison == "swhu_vs_huhu" ~ "swine-human > human-human",
      comparison == "huhu_vs_swhu" ~ "swine-human < human-human",
      comparison == "husw_vs_huhu" ~ "human-swine > human-human",
      comparison == "huhu_vs_husw" ~ "human-swine < human-human",
      comparison == "swhu_vs_swsw" ~ "swine-human > swine-swine",
      comparison == "swsw_vs_swhu" ~ "swine-human < swine-swine")
  )

within_host_data <- summary %>%
  filter(comparison %in% within_host_comparisons) %>%
  mutate(
    comparison = case_when(
      comparison == "huhu_vs_swsw" ~ "human-human > swine-swine",
      comparison == "swsw_vs_huhu" ~ "human-human < swine-swine")
  )

comparison_order <- c("swine-human > human-human", "swine-human < human-human",
                      "human-swine > human-human","human-swine < human-human",
                      "swine-human > swine-swine", "swine-human < swine-swine",
                      "human-swine > swine-swine", "human-swine < swine-swine",
                      "human-human > swine-swine", "human-human < swine-swine")

cross_host_data$comparison <- factor(cross_host_data$comparison, levels=comparison_order)

within_host_data$comparison <- factor(within_host_data$comparison, levels=comparison_order)


comparison_colors <- c(
  "human-swine > swine-swine" = "#89B177",
  "human-swine < swine-swine" = "#89B177",
  "human-swine > human-human" = "#C1E3AE",
  "human-swine < human-human" = "#C1E3AE",
  "swine-human > swine-swine" = "#FFE6AA",
  "swine-human < swine-swine" = "#FFE6AA",
  "swine-human > human-human" = "#FFBF60",
  "swine-human < human-human" = "#FFBF60",
  "human-human > swine-swine" = "#D78B5E",
  "human-human < swine-swine" = "#45BACF")

plot_c <- ggplot(cross_host_data, aes(x = comparison, y = count, fill = comparison, linetype = comparison)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9), 
           color = "black", 
           linewidth = 0.25) +
  labs(title = "",
       x = "",
       y = "Number of positions with higher substitution rate") +
  scale_fill_manual(values = comparison_colors) +
  scale_linetype_manual(values = c(
    "human-swine > swine-swine" = "solid",
    "human-swine < swine-swine" = "dashed",
    "swine-human > human-human" = "solid",
    "swine-human < human-human" = "dashed",
    "human-swine > human-human" = "solid",
    "human-swine < human-human" = "dashed",
    "swine-human > swine-swine" = "solid",
    "swine-human < swine-swine" = "dashed"
  )) +
  facet_wrap(~segment, ncol = 7) +
  theme_classic() +
  theme1
plot_c <- plot_c + guides(fill = guide_legend(nrow = 4, byrow = TRUE))


plot_d <- ggplot(within_host_data, aes(x = segment, y = count, fill = comparison, linetype = comparison)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.85), 
           color = "black",
           linewidth = 0.25) +
  labs(title = "",
       x = "",
       y = "Number of positions with higher substitution rate") +
  scale_fill_manual(values = comparison_colors) +
  scale_linetype_manual(values = c(
    "human-human > swine-swine" = "solid",
    "human-human < swine-swine" = "solid")) +
  theme_classic() +
  theme3


combined_plot <- ggarrange(plot_c, plot_d,
                           ncol = 1,
                           nrow = 2,
                           heights = c(2,1),
                           labels = c("c", "d"),
                           font.label = list(size = 10))


combined_plot_final <- ggarrange(total_and_normalized_mut_count_plots, combined_plot,
                           ncol = 2,
                           nrow = 1)

ggsave("figures_tables/figure2_mut_count_plots/Figure2_mut_count_and_rate_plots.png", plot = combined_plot_final, dpi = 600, width = 18 , height = 20, units = "cm" )


###############################################################################


bayes_positions <- combined_data %>%
  filter(comparison %in% c("swhu_vs_huhu", "husw_vs_swsw")) %>%
  select(segment, seqpos, comparison) %>%
  dplyr::rename(position = seqpos, trait = comparison) %>%
  mutate(trait = recode(trait, 
                        "swhu_vs_huhu" = "swine-human", 
                        "husw_vs_swsw" = "human-swine"))

write.table(bayes_positions, "figures_tables/figure2_mut_count_plots/bayes_positions.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

