library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotrix)

setwd("PATH_TO_PROJECT_DIRECTORY")

###############################################################
# Load data
pb2_fubar <- read.csv("selection_analysis/pb2/pb2_fubar.csv")
colnames(pb2_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

pb2_fubar1 <- pb2_fubar %>%
  mutate(omega = beta/alpha,
         protein = "PB2")

pb1_fubar <- read.csv("selection_analysis/pb1/pb1_fubar.csv")
colnames(pb1_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

pb1_fubar1 <- pb1_fubar %>%
  mutate(omega = beta/alpha,
         protein = "PB1")

pb1_f2_fubar <- read.csv("selection_analysis/pb1/pb1-f2_fubar.csv")
colnames(pb1_f2_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

pb1_f2_fubar1 <- pb1_f2_fubar %>%
  mutate(omega = beta/alpha,
         protein = "PB1-F2")

pa_fubar <- read.csv("selection_analysis/pa/pa_fubar.csv")
colnames(pa_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

pa_fubar1 <- pa_fubar %>%
  mutate(omega = beta/alpha,
         protein = "PA")

pa_x_fubar <- read.csv("selection_analysis/pa/pa-x_fubar.csv")
colnames(pa_x_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

pa_x_fubar1 <- pa_x_fubar %>%
  mutate(omega = beta/alpha,
         protein = "PA-X")

h1_fubar <- read.csv("selection_analysis/ha/h1_fubar.csv")
colnames(h1_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

h1_fubar1 <- h1_fubar %>%
  mutate(omega = beta/alpha,
         protein = "H1")

h3_fubar <- read.csv("selection_analysis/ha/h3_fubar.csv")
colnames(h3_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

h3_fubar1 <- h3_fubar %>%
  mutate(omega = beta/alpha,
         protein = "H3")

np_fubar <- read.csv("selection_analysis/np/np_fubar.csv")
colnames(np_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

np_fubar1 <- np_fubar %>%
  mutate(omega = beta/alpha,
         protein = "NP")

n1_fubar <- read.csv("selection_analysis/na/n1_fubar.csv")
colnames(n1_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

n1_fubar1 <- n1_fubar %>%
  mutate(omega = beta/alpha,
         protein = "N1")

n2_fubar <- read.csv("selection_analysis/na/n2_fubar.csv")
colnames(n2_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

n2_fubar1 <- n2_fubar %>%
  mutate(omega = beta/alpha,
         protein = "N2")

m1_fubar <- read.csv("selection_analysis/mp/m1_fubar.csv")
colnames(m1_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

m1_fubar1 <- m1_fubar %>%
  mutate(omega = beta/alpha,
         protein = "M1")

m2_fubar <- read.csv("selection_analysis/mp/m2_fubar.csv")
colnames(m2_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

m2_fubar1 <- m2_fubar %>%
  mutate(omega = beta/alpha,
         protein = "M2")

ns1_fubar <- read.csv("selection_analysis/ns/ns1_fubar.csv")
colnames(ns1_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

ns1_fubar1 <- ns1_fubar %>%
  mutate(omega = beta/alpha,
         protein = "NS1")

nep_fubar <- read.csv("selection_analysis/ns/nep_fubar.csv")
colnames(nep_fubar) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha", "prob_neg", "prob_pos", "BayesFactor[alpha<beta]")

nep_fubar1 <- nep_fubar %>%
  mutate(omega = beta/alpha,
         protein = "NEP")

combined <- rbind(pb2_fubar1, pb1_fubar1, pb1_f2_fubar1, pa_fubar1, pa_x_fubar1,
                  h1_fubar1, h3_fubar1, np_fubar1, n1_fubar1, n2_fubar1, 
                  m1_fubar1, m2_fubar1, ns1_fubar1, nep_fubar1)

combined$protein <- factor(combined$protein, levels=c("PB2", "PB1", "PB1-F2", "PA",
                                                      "PA-X", "H1", "H3", "NP",
                                                      "N1", "N2", "M1", "M2",
                                                      "NS1", "NEP"))

combined_summarized <- combined %>%
  group_by(protein) %>%
  summarize(average_omega = mean(omega, na.rm = TRUE), 
            n_neg = sum(prob_neg>0.9), 
            n_pos = sum(prob_pos>0.9),
            total = n()) %>%
  mutate("purifying" = (n_neg/total)*100,
         "positive" = (n_pos/total)*100)


combined_long <- combined_summarized %>%
  pivot_longer(cols = c(purifying, positive), names_to = "type", values_to = "proportion")

# Plotting

theme1 <- theme(legend.position = "top",
                legend.title = element_blank(),
                legend.key.size = unit(7, "pt"),
                legend.text = element_text(size = 7),
                legend.spacing = unit(-0.5, "pt"),
                legend.margin = margin(unit(c(0,0,10,0), "lines")),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 7),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 7),
                axis.ticks.y = element_line(size = 0.2),
                axis.ticks.x = element_blank(),
                strip.background = element_rect(linewidth = 0.5),
                strip.text = element_text(size = 8),
                axis.line = element_line(size = 0.1),
                panel.spacing.x = unit(0.1, "lines"), 
                panel.spacing.y = unit(0.2, "lines"),
                panel.border = element_rect(fill = NA, color = "grey10", linewidth = 0.5),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                legend.box.spacing = unit(0.1, "lines"),
                plot.margin = unit(c(0.2,0.5,1,0.5), "lines"))

theme2 <- theme(legend.position = "top",
                legend.title = element_blank(),
                legend.key.size = unit(7, "pt"),
                legend.text = element_text(size = 7),
                legend.spacing = unit(-0.5, "pt"),
                legend.margin = margin(unit(c(0,0,10,0), "lines")),
                axis.text.x = element_text(size = 7, vjust = 0.8, angle = 45, hjust = 0.9),
                axis.text.y = element_text(size=7),
                axis.title.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.ticks = element_line(size = 0.2),
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

proportions <- ggplot(combined_long) +
  geom_col(aes(x = type, y = proportion, fill = type), position = "dodge") +
  labs(title = "",
       x = "",
       y = "% codons under selection") +
  scale_fill_manual(values = c("#3B5A9B","#A43A3A")) +
  facet_wrap(~protein, nrow=1) +
  theme_classic() +
  theme1

distributions <- ggplot(combined) +
  geom_point(aes(x = site, y = omega), color = "grey20", size = 0.5) +
  geom_line(aes(x = site, y = omega), color = "grey20", size = 0.2) +
  geom_point(data = subset(combined, prob_neg > 0.9), aes(x = site, y = omega), color = "#A43A3A", size = 0.5) +
  geom_point(data = subset(combined, prob_pos > 0.9), aes(x = site, y = omega), color = "#3B5A9B", size = 0.5) +
  labs(title = "",
       x = "position",
       y = "ω (dN/dS)") +
  facet_wrap(~protein, ncol = 4, scales = "free_x") +
  scale_x_continuous(n.breaks = 10) +
  theme_classic() +
  theme2



FigureS3_combined_plot <- ggarrange(proportions, distributions,
                                   ncol = 1,
                                   nrow = 2,
                                   heights = c(1,3),
                                   common.legend = TRUE,
                                   labels = c("a", "b"), 
                                   font.label = list(size = 10))

ggsave("figures_tables/figure3_S3S4_selection/figureS3_fubar_overview.png", 
       plot = FigureS3_combined_plot, dpi = 300, width = 17 , height = 25, units = "cm", bg = "white")




