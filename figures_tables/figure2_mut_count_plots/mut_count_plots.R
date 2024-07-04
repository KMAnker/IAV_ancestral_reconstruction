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
                legend.key.size = unit(9, "pt"),
                legend.text = element_text(size = 8),
                legend.spacing = unit(-0.5, "pt"),
                legend.margin = margin(unit(c(0.5,0,5,0), "lines")),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=8),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 8),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                strip.background = element_rect(linewidth = 0.5),
                strip.text = element_text(size = 9),
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

theme2 <- theme(legend.position = "right",
                legend.title = element_blank(),
                legend.key.size = unit(8, "pt"),
                legend.text = element_text(size = 8),
                legend.spacing = unit(-0.5, "pt"),
                legend.margin = margin(unit(c(0.1,0,0.5,35), "lines")),
                axis.text.x = element_text(size = 9, vjust = 0.8, angle = 45, hjust = 0.9),
                axis.text.y = element_text(size=8),
                axis.title.x = element_text(size = 8),
                axis.title.y = element_text(size = 8),
                axis.ticks.y = element_blank(),
                strip.background = element_rect(linewidth = 0.5),
                strip.text = element_text(size = 9),
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
ggsave("figures_tables/figure2_mut_count_plots/Figure2_total_and_normalized_mut_count_plots.png", plot = total_and_normalized_mut_count_plots, dpi = 300, width = 12 , height = 18, units = "cm" )

