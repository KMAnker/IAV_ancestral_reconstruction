library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)

setwd("PATH/TO/PROJECT/DIRECTORY")

auc_table <- read.table("AUC_table.txt", header = TRUE, na.strings = "")

# Reshape data from wide to long format
df_long <- auc_table %>%
  pivot_longer(
    cols = starts_with("AUC_"),
    names_to = "AUC_type",
    values_to = "AUC_value"
  ) %>%
  mutate(
    AUC_type = recode(AUC_type,
                      "AUC_hu.hu" = "human-human",
                      "AUC_sw.sw" = "swine-swine",
                      "AUC_sw.hu" = "swine-human",
                      "AUC_hu.sw" = "human-swine"))

df_long$protein <- factor(df_long$protein, levels= c("PB2", "PB1", "PB1-F2","PA", "PA-X", "H1", "H3", "NP", "N1", "N2", "M1", "M2", "NS1", "NEP"))
df_long$AUC_type <- factor(df_long$AUC_type, levels= c("human-human", "swine-swine", "human-swine", "swine-human"))

# Make a heatmap using ggplot
theme1 <- theme(axis.text.x = element_text(size = 8),
                axis.text.y = element_text(size = 8),
                panel.grid = element_blank(),
                legend.box.spacing = unit(0.1, "lines"),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 7), 
                axis.title.x = element_text(size = 8, vjust = -2),
                axis.title.y = element_text(size = 8),
                plot.margin = margin(5, 5, 10, 5), "points")


plot <- ggplot(df_long, aes(x = AUC_type, y = protein, fill = AUC_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", AUC_value)), size = 3) +
  scale_y_discrete(limits=rev) +
  scale_fill_distiller(palette = "RdBu", name = "AUC", limits = c(0, 1), direction = 1) +
  theme_minimal() +
  labs(x = "Transmission type", y = "Protein", fill = "AUC") +
  theme1

ggsave("figures_tables/figure4_5_5_S4_S5_xgboost/AUC_heatmap.png", plot = plot, dpi = 600, width = 12 , height = 12, units = "cm", bg = "white")
