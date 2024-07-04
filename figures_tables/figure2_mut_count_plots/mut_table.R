library(tidyverse)
library(dplyr)

setwd("PATH_TO_PROJECT_DIRECTORY")


## loading data

pb2_aa_branchdiff <- read.table("anclib_files/pb2/pb2_branchdiffinfo_aa_br2.txt")
colnames(pb2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pb1_aa_branchdiff <- read.table("anclib_files/pb1/pb1_branchdiffinfo_aa_br2.txt")
colnames(pb1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pb1_f2_aa_branchdiff <- read.table("anclib_files/pb1/pb1-f2_branchdiffinfo_aa_br2.txt")
colnames(pb1_f2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pa_aa_branchdiff <- read.table("anclib_files/pa/pa_branchdiffinfo_aa_br2.txt")
colnames(pa_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

pa_x_aa_branchdiff <- read.table("anclib_files/pa/pa-x_branchdiffinfo_aa_br2.txt")
colnames(pa_x_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

h1_aa_branchdiff <- read.table("anclib_files/ha/h1_branchdiffinfo_aa_br2.txt")
colnames(h1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

h3_aa_branchdiff <- read.table("anclib_files/ha/h3_branchdiffinfo_aa_br2.txt")
colnames(h3_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

np_aa_branchdiff <- read.table("anclib_files/np/np_branchdiffinfo_aa_br2.txt")
colnames(np_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

n1_aa_branchdiff <- read.table("anclib_files/na/n1_branchdiffinfo_aa_br2.txt")
colnames(n1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

n2_aa_branchdiff <- read.table("anclib_files/na/n2_branchdiffinfo_aa_br2.txt")
colnames(n2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

m1_aa_branchdiff <- read.table("anclib_files/mp/m1_branchdiffinfo_aa_br2.txt")
colnames(m1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

m2_aa_branchdiff <- read.table("anclib_files/mp/m2_branchdiffinfo_aa_br2.txt")
colnames(m2_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

ns1_aa_branchdiff <- read.table("anclib_files/ns/ns1_branchdiffinfo_aa_br2.txt")
colnames(ns1_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")

nep_aa_branchdiff <- read.table("anclib_files/ns/nep_branchdiffinfo_aa_br2.txt")
colnames(nep_aa_branchdiff) <- c("node_from", "node_to", "branchlen", "seq_pos", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "branch_type", "external")



# Function to process the branchdiff data frames and create the summarized mutation data
process_branchdiff <- function(data, protein_name) {
  data %>%
    filter(branchlen <= 15) %>%
    mutate(trait = paste(trait_from, trait_to, sep = "-"), .before = node_from,
           protein = protein_name) %>%
    filter(residue_from != 'X', residue_to != 'X', residue_from != '*', residue_to != '*') %>%
    group_by(trait, seq_pos) %>%
    filter(residue_from != residue_to) %>%
    mutate(mut = paste(residue_from, residue_to, sep = "-")) %>%
    ungroup() %>%
    select(protein, trait, seq_pos, mut) %>%
    group_by(mut, trait) %>%
    mutate(mut_count = n()) %>%
    ungroup()
}

# Function to calculate trait counts
calculate_trait_counts <- function(data, protein_name) {
  data %>%
    filter(branchlen <= 15) %>%
    mutate(trait = paste(trait_from, trait_to, sep = "-"), .before = node_from,
           branch = paste(node_from, node_to, sep = "-"),
           protein = protein_name) %>%
    group_by(protein, trait) %>%
    summarize(num_samples = n_distinct(branch)) %>%
    mutate(total_num_samples = sum(num_samples)) %>%
    ungroup()
}

# Function to combine the results
combine_results <- function(branchdiff_data, protein_name) {
  summarised_data <- process_branchdiff(branchdiff_data, protein_name)
  trait_counts <- calculate_trait_counts(branchdiff_data, protein_name)
  
  summarised_data %>%
    group_by(protein, trait, seq_pos, mut) %>%
    summarise(mutfreq = n()) %>%
    ungroup() %>%
    left_join(trait_counts, by = c("protein", "trait")) %>%
    mutate(rel_mutfreq = (mutfreq / num_samples) * 100)
}


# Combine all branchdiff files to a list
combined_branchdiff_files <- list(PB2 = pb2_aa_branchdiff, PB1 = pb1_aa_branchdiff, "PB1-F2" = pb1_f2_aa_branchdiff,
                                  PA = pa_aa_branchdiff, "PA-X" = pa_x_aa_branchdiff,  H1 = h1_aa_branchdiff,
                                  H3 = h3_aa_branchdiff, NP = np_aa_branchdiff, N1 = n1_aa_branchdiff, 
                                  N2 = n2_aa_branchdiff, M1 = m1_aa_branchdiff, M2 = m2_aa_branchdiff,
                                  NS1 = ns1_aa_branchdiff, NEP = nep_aa_branchdiff)

# Apply the function to each data frame in the list
summarised_mut_results <- lapply(names(combined_branchdiff_files), function(protein_name) {
  combine_results(combined_branchdiff_files[[protein_name]], protein_name)
})


# Combine all results into a single data frame
mut_table <- bind_rows(summarised_mut_results)


write.table(mut_table, "figures_tables/figure2_mut_count_plots/mut_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
