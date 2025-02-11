library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)

setwd("PATH/TO/PROJECT/DIRECTORY")

# Function to reshape dataframes
process_shap_importance <- function(shap_importance_df) {
  # Gather the dataframe into key-value pairs
  long_df <- shap_importance_df %>%
    tidyr::gather(key = "key", value = "value", -Feature)
  
  # Process the dataframe to extract traits, types, and presence
  long_df <- long_df %>%
    dplyr::mutate(
      trait = stringr::str_extract(key, "^(human_human|swine_swine|human_swine|swine_human)"),
      type = dplyr::case_when(
        stringr::str_detect(key, "_HOW_important$") ~ "HOW_important",
        stringr::str_detect(key, "_importance_not_present$") ~ "importance_not_present",
        stringr::str_detect(key, "_std_not_present$") ~ "std_not_present",
        stringr::str_detect(key, "_not_present_N$") ~ "not_present_N",
        stringr::str_detect(key, "_importance_present$") ~ "importance_present",
        stringr::str_detect(key, "_std_present$") ~ "std_present",
        stringr::str_detect(key, "_present_N$") ~ "present_N"
      ),
      presence = dplyr::case_when(
        stringr::str_detect(key, "_not_present_N$") ~ "not_present",
        stringr::str_detect(key, "_not_present$") ~ "not_present",
        stringr::str_detect(key, "_present_N$") ~ "present",
        stringr::str_detect(key, "_present$") ~ "present",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-key)
  
  # Split the dataframe into separate dataframes based on type
  how_important_df <- long_df %>% dplyr::filter(type == "HOW_important") %>% dplyr::select(Feature, trait, HOW_important = value)
  importance_df <- long_df %>% dplyr::filter(stringr::str_detect(type, "importance")) %>% dplyr::select(Feature, trait, presence, importance = value)
  std_df <- long_df %>% dplyr::filter(stringr::str_detect(type, "std")) %>% dplyr::select(Feature, trait, presence, std = value)
  n_df <- long_df %>% dplyr::filter(stringr::str_detect(type, "_N$")) %>% dplyr::select(Feature, trait, presence, N = value)
  
  # Join the dataframes together
  final_df <- importance_df %>%
    dplyr::left_join(std_df, by = c("Feature", "trait", "presence")) %>%
    dplyr::left_join(n_df, by = c("Feature", "trait", "presence"))
}

normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

traits <- c("human_human", "swine_swine", "human_swine", "swine_human")

trait_order <- c("human_human","swine_swine","human_swine", "swine_human")


############

# PB2
shap_importance_pb2 <- read.csv("xgboost/pb2/pb2_importance_dataframe.csv", header = TRUE)

pb2_plotting_df <- process_shap_importance(shap_importance_pb2)

pb2_plotting_df$trait <- factor(pb2_plotting_df$trait, levels=trait_order)
pb2_plotting_df$presence <- factor(pb2_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_pb2 <- unique(pb2_plotting_df$Feature)
ordered_levels_pb2 <- ordered_levels_pb2[order(as.numeric(str_extract(ordered_levels_pb2, "\\d+")))]
pb2_plotting_df$Feature <- factor(pb2_plotting_df$Feature, levels = ordered_levels_pb2)


top5_pb2 <- pb2_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


pb2 <- ggplot(top5_pb2, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "PB2",
       x = "",
       y = "SHAP value")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        legend.spacing = unit(0, "lines"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


pb2_positions <- top5_pb2 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "PB2") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

pb2_positions_complete <- expand.grid(
  position = unique(pb2_positions$position),
  trait = traits,
  segment = "PB2"
)

# PB1
shap_importance_pb1 <- read.csv("xgboost/pb1/pb1_importance_dataframe.csv", header = TRUE)

pb1_plotting_df <- process_shap_importance(shap_importance_pb1)

pb1_plotting_df$trait <- factor(pb1_plotting_df$trait, levels=trait_order)
pb1_plotting_df$presence <- factor(pb1_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_pb1 <- unique(pb1_plotting_df$Feature)
ordered_levels_pb1 <- ordered_levels_pb1[order(as.numeric(str_extract(ordered_levels_pb1, "\\d+")))]
pb1_plotting_df$Feature <- factor(pb1_plotting_df$Feature, levels = ordered_levels_pb1)


top5_pb1 <- pb1_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


pb1 <- ggplot(top5_pb1, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "PB1",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))
  


pb1_positions <- top5_pb1 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "PB1") %>%
  select(segment, position, trait) %>%
  distinct()%>%
  arrange(position)

pb1_positions_complete <- expand.grid(
  position = unique(pb1_positions$position),
  trait = traits,
  segment = "PB1"
)

# PB1-F2
shap_importance_pb1f2 <- read.csv("xgboost/pb1/pb1-f2_importance_dataframe.csv", header = TRUE)

pb1f2_plotting_df <- process_shap_importance(shap_importance_pb1f2)

pb1f2_plotting_df$trait <- factor(pb1f2_plotting_df$trait, levels=trait_order)
pb1f2_plotting_df$presence <- factor(pb1f2_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_pb1f2 <- unique(pb1f2_plotting_df$Feature)
ordered_levels_pb1f2 <- ordered_levels_pb1f2[order(as.numeric(str_extract(ordered_levels_pb1f2, "\\d+")))]
pb1f2_plotting_df$Feature <- factor(pb1f2_plotting_df$Feature, levels = ordered_levels_pb1f2)


top5_pb1f2 <- pb1f2_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


pb1f2 <- ggplot(top5_pb1f2, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "PB1-F2",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))
  

pb1f2_positions <- top5_pb1f2 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "PB1-F2") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

pb1f2_positions_complete <- expand.grid(
  position = unique(pb1f2_positions$position),
  trait = traits,
  segment = "PB1-F2"
)

# PA
shap_importance_pa <- read.csv("xgboost/pa/pa_importance_dataframe.csv", header = TRUE)

pa_plotting_df <- process_shap_importance(shap_importance_pa)

pa_plotting_df$trait <- factor(pa_plotting_df$trait, levels=trait_order)
pa_plotting_df$presence <- factor(pa_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_pa <- unique(pa_plotting_df$Feature)
ordered_levels_pa <- ordered_levels_pa[order(as.numeric(str_extract(ordered_levels_pa, "\\d+")))]
pa_plotting_df$Feature <- factor(pa_plotting_df$Feature, levels = ordered_levels_pa)


top5_pa <- pa_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


pa <- ggplot(top5_pa, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "PA",
       x = "",
       y = "SHAP value")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


pa_positions <- top5_pa %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "PA") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

pa_positions_complete <- expand.grid(
  position = unique(pa_positions$position),
  trait = traits,
  segment = "PA"
)

# PA-X
shap_importance_pax <- read.csv("xgboost/pa/pa-x_importance_dataframe.csv", header = TRUE)

pax_plotting_df <- process_shap_importance(shap_importance_pax)

pax_plotting_df$trait <- factor(pax_plotting_df$trait, levels=trait_order)
pax_plotting_df$presence <- factor(pax_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_pax <- unique(pax_plotting_df$Feature)
ordered_levels_pax <- ordered_levels_pax[order(as.numeric(str_extract(ordered_levels_pax, "\\d+")))]
pax_plotting_df$Feature <- factor(pax_plotting_df$Feature, levels = ordered_levels_pax)


top5_pax <- pax_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


pax <- ggplot(top5_pax, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "PA-X",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


pax_positions <- top5_pax %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "PA-X") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

pax_positions_complete <- expand.grid(
  position = unique(pax_positions$position),
  trait = traits,
  segment = "PA-X"
)

# H1
shap_importance_h1 <- read.csv("xgboost/ha/h1_importance_dataframe.csv", header = TRUE)

h1_plotting_df <- process_shap_importance(shap_importance_h1)

h1_plotting_df$trait <- factor(h1_plotting_df$trait, levels=trait_order)
h1_plotting_df$presence <- factor(h1_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_h1 <- unique(h1_plotting_df$Feature)
ordered_levels_h1 <- ordered_levels_h1[order(as.numeric(str_extract(ordered_levels_h1, "\\d+")))]
h1_plotting_df$Feature <- factor(h1_plotting_df$Feature, levels = ordered_levels_h1)


top5_h1 <- h1_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)

  
h1 <- ggplot(top5_h1, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "H1",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))
  

h1_positions <- top5_h1 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "H1") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

h1_positions_complete <- expand.grid(
  position = unique(h1_positions$position),
  trait = traits,
  segment = "H1"
)


# H3
shap_importance_h3 <- read.csv("xgboost/ha/h3_importance_dataframe.csv", header = TRUE)

h3_plotting_df <- process_shap_importance(shap_importance_h3)

h3_plotting_df$trait <- factor(h3_plotting_df$trait, levels=trait_order)
h3_plotting_df$presence <- factor(h3_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_h3 <- unique(h3_plotting_df$Feature)
ordered_levels_h3 <- ordered_levels_h3[order(as.numeric(str_extract(ordered_levels_h3, "\\d+")))]
h3_plotting_df$Feature <- factor(h3_plotting_df$Feature, levels = ordered_levels_h3)


top5_h3 <- h3_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


h3 <- ggplot(top5_h3, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "H3",
       x = "",
       y = "SHAP value")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines")) 
  

h3_positions <- top5_h3 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "H3") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

h3_positions_complete <- expand.grid(
  position = unique(h3_positions$position),
  trait = traits,
  segment = "H3"
)

# NP
shap_importance_np <- read.csv("xgboost/np/np_importance_dataframe.csv", header = TRUE)

np_plotting_df <- process_shap_importance(shap_importance_np)

np_plotting_df$trait <- factor(np_plotting_df$trait, levels=trait_order)
np_plotting_df$presence <- factor(np_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_np <- unique(np_plotting_df$Feature)
ordered_levels_np <- ordered_levels_np[order(as.numeric(str_extract(ordered_levels_np, "\\d+")))]
np_plotting_df$Feature <- factor(np_plotting_df$Feature, levels = ordered_levels_np)


top5_np <- np_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)

  
np <- ggplot(top5_np, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "NP",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))  


np_positions <- top5_np %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "NP") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

np_positions_complete <- expand.grid(
  position = unique(np_positions$position),
  trait = traits,
  segment = "NP"
)

# N1
shap_importance_n1 <- read.csv("xgboost/na/n1_importance_dataframe.csv", header = TRUE)

n1_plotting_df <- process_shap_importance(shap_importance_n1)

n1_plotting_df$trait <- factor(n1_plotting_df$trait, levels=trait_order)
n1_plotting_df$presence <- factor(n1_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_n1 <- unique(n1_plotting_df$Feature)
ordered_levels_n1 <- ordered_levels_n1[order(as.numeric(str_extract(ordered_levels_n1, "\\d+")))]
n1_plotting_df$Feature <- factor(n1_plotting_df$Feature, levels = ordered_levels_n1)


top5_n1 <- n1_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


n1 <- ggplot(top5_n1, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "N1",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines")) 


n1_positions <- top5_n1 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "N1") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

n1_positions_complete <- expand.grid(
  position = unique(n1_positions$position),
  trait = traits,
  segment = "N1"
)

# N2
shap_importance_n2 <- read.csv("xgboost/na/n2_importance_dataframe.csv", header = TRUE)

n2_plotting_df <- process_shap_importance(shap_importance_n2)


n2_plotting_df$trait <- factor(n2_plotting_df$trait, levels=trait_order)
n2_plotting_df$presence <- factor(n2_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_n2 <- unique(n2_plotting_df$Feature)
ordered_levels_n2 <- ordered_levels_n2[order(as.numeric(str_extract(ordered_levels_n2, "\\d+")))]
n2_plotting_df$Feature <- factor(n2_plotting_df$Feature, levels = ordered_levels_n2)


top5_n2 <- n2_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


  n2 <- ggplot(top5_n2, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "N2",
       x = "",
       y = "SHAP value")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))
  

n2_positions <- top5_n2 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "N2") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

n2_positions_complete <- expand.grid(
  position = unique(n2_positions$position),
  trait = traits,
  segment = "N2"
)

# M1
shap_importance_m1 <- read.csv("xgboost/mp/m1_importance_dataframe.csv", header = TRUE)

m1_plotting_df <- process_shap_importance(shap_importance_m1)

m1_plotting_df$trait <- factor(m1_plotting_df$trait, levels=trait_order)
m1_plotting_df$presence <- factor(m1_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_m1 <- unique(m1_plotting_df$Feature)
ordered_levels_m1 <- ordered_levels_m1[order(as.numeric(str_extract(ordered_levels_m1, "\\d+")))]
m1_plotting_df$Feature <- factor(m1_plotting_df$Feature, levels = ordered_levels_m1)


top5_m1 <- m1_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


m1 <- ggplot(top5_m1, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "M1",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


m1_positions <- top5_m1 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "M1") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

m1_positions_complete <- expand.grid(
  position = unique(m1_positions$position),
  trait = traits,
  segment = "M1"
)

# M2
shap_importance_m2 <- read.csv("xgboost/mp/m2_importance_dataframe.csv", header = TRUE)

m2_plotting_df <- process_shap_importance(shap_importance_m2)

m2_plotting_df$trait <- factor(m2_plotting_df$trait, levels=trait_order)
m2_plotting_df$presence <- factor(m2_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_m2 <- unique(m2_plotting_df$Feature)
ordered_levels_m2 <- ordered_levels_m2[order(as.numeric(str_extract(ordered_levels_m2, "\\d+")))]
m2_plotting_df$Feature <- factor(m2_plotting_df$Feature, levels = ordered_levels_m2)


top5_m2 <- m2_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


m2 <- ggplot(top5_m2, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "M2",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


m2_positions <- top5_m2 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "M2") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

m2_positions_complete <- expand.grid(
  position = unique(m2_positions$position),
  trait = traits,
  segment = "M2"
)

# NS1
shap_importance_ns1 <- read.csv("xgboost/ns/ns1_importance_dataframe.csv", header = TRUE)

ns1_plotting_df <- process_shap_importance(shap_importance_ns1)

ns1_plotting_df$trait <- factor(ns1_plotting_df$trait, levels=trait_order)
ns1_plotting_df$presence <- factor(ns1_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_ns1 <- unique(ns1_plotting_df$Feature)
ordered_levels_ns1 <- ordered_levels_ns1[order(as.numeric(str_extract(ordered_levels_ns1, "\\d+")))]
ns1_plotting_df$Feature <- factor(ns1_plotting_df$Feature, levels = ordered_levels_ns1)


top5_ns1 <- ns1_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)



ns1 <- ggplot(top5_ns1, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "NS1",
       x = "",
       y = "SHAP value")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


ns1_positions <- top5_ns1 %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "NS1") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

ns1_positions_complete <- expand.grid(
  position = unique(ns1_positions$position),
  trait = traits,
  segment = "NS1"
)

# NEP
shap_importance_nep <- read.csv("xgboost/ns/nep_importance_dataframe.csv", header = TRUE)

nep_plotting_df <- process_shap_importance(shap_importance_nep)

nep_plotting_df$trait <- factor(nep_plotting_df$trait, levels=trait_order)
nep_plotting_df$presence <- factor(nep_plotting_df$presence, levels=c("present", "not_present"))

ordered_levels_nep <- unique(nep_plotting_df$Feature)
ordered_levels_nep <- ordered_levels_nep[order(as.numeric(str_extract(ordered_levels_nep, "\\d+")))]
nep_plotting_df$Feature <- factor(nep_plotting_df$Feature, levels = ordered_levels_nep)


top5_nep <- nep_plotting_df %>%
  mutate(norm_N = normalize(N)) %>%
  group_by(trait) %>%
  arrange(desc(abs(importance))) %>%
  slice_head(n = 5)


nep <- ggplot(top5_nep, aes(x = Feature, y = importance, fill = trait, alpha = norm_N)) +
  geom_bar(aes(fill = trait, linetype = presence), stat = "identity", 
           position = "dodge", width = 0.9, color = "grey40", size = 0.3) +
  scale_fill_manual(values = c("#D78B5E","#45BACF", "#A5CD92", "#FFD685"),
                    labels= c(human_human="human-human", swine_swine="swine-swine", human_swine="human-swine", swine_human="swine-human")) +
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
  labs(title = "NEP",
       x = "",
       y = "")  +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(7, "points"),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.9, hjust = 0.9),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6),
        panel.spacing.x = unit(0.2, "lines"), 
        panel.spacing.y = unit(0.1, "lines"),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"))


nep_positions <- top5_nep %>%
  mutate(position = as.numeric(str_extract(Feature, "\\d+")),
         segment = "NEP") %>%
  select(segment, position, trait) %>%
  distinct() %>%
  arrange(position)

nep_positions_complete <- expand.grid(
  position = unique(nep_positions$position),
  trait = traits,
  segment = "NEP"
)

##########

plot_all <- ggarrange(pb2, pb1, pb1f2, pa, pax, h1, h3, np, n1, n2, m1, m2, ns1, nep, 
          ncol = 3, nrow = 5,
          common.legend = TRUE, 
          legend = "top",
          font.label = list(size = 10))

ggsave("figures_tables/figure4_5_S4_S5_xgboost/Figure5_shap_values_plot.png", plot = plot_all, dpi = 600, width = 17 , height = 25, units = "cm", bg = "white")


####################################################

positions_all <- rbind(pb2_positions_complete, pb1_positions_complete, pb1f2_positions_complete, pa_positions_complete, pax_positions_complete, 
                       h1_positions_complete, h3_positions_complete, np_positions_complete, n1_positions_complete, n2_positions_complete, 
                       m1_positions_complete, m2_positions_complete, ns1_positions_complete, nep_positions_complete) %>%
  mutate(trait = as.character(trait),
         trait = str_replace(trait, "_", "-"))%>%
  ungroup()

positions_all_2 <- rbind(pb2_positions, pb1_positions, pb1f2_positions, pa_positions, pax_positions, 
                       h1_positions, h3_positions, np_positions, n1_positions, n2_positions, 
                       m1_positions, m2_positions, ns1_positions, nep_positions) %>%
  mutate(trait = as.character(trait),
         trait = str_replace(trait, "_", "-"))%>%
  ungroup()

write.table(positions_all_2, "figures_tables/figure4_S4_S5_xgboost/shap_positions.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################

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

summarize_branch_diff <- function(file_path, segment_name) {
  # Read the data
  branchdiff <- read.table(file_path)
  colnames(branchdiff) <- c("node_from", "node_to", "branchlen", "position", "trait_from", "trait_to", "traitprob_from", "traitprob_to", "residue_from", "residue_to", "trait", "external")
  
  # Calculate total residue counts for all rows
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
    arrange(desc(n_residue), .by_group = TRUE) %>%
    summarise(
      "Amino Acids (%)" = paste0(residue_to, "(", round(n_residue/total_residue * 100, 1), "%)", collapse = ", ")
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
    arrange(desc(n_mutation), .by_group = TRUE) %>%
    summarise(
      "Mutations (%)" = ifelse(sum(n_mutation) == 0, "NA",
                               paste0(mutation, "(", round(n_mutation/total_residue * 100, 1), "%)", collapse = ", "))
    ) %>%
    ungroup() %>%
    mutate(segment = segment_name, .before = "position")
  
  # Calculate total residue counts for human and swine observed cases
  total_human_observed_counts <- branchdiff %>%
    filter(external == 1, trait_to == "human", branchlen <= 15) %>%
    mutate(segment = segment_name,
           position = as.numeric(position)) %>%
    group_by(position) %>%
    summarise(total_human_observed = n()) %>%
    ungroup()
  
  total_swine_observed_counts <- branchdiff %>%
    filter(external == 1, trait_to == "swine", branchlen <= 15) %>%
    mutate(segment = segment_name,
           position = as.numeric(position)) %>%
    group_by(position) %>%
    summarise(total_swine_observed = n()) %>%
    ungroup()
  
  # Tally for human and swine observed AAs
  human_observed_residue_tally <- branchdiff %>%
    filter(external == 1, trait_to == "human", branchlen <= 15) %>%
    mutate(position = as.numeric(position)) %>%
    group_by(position, residue_to) %>%
    tally() %>%
    dplyr::rename(n_human_residue = n) %>%
    ungroup()
  
  swine_observed_residue_tally <- branchdiff %>%
    filter(external == 1, trait_to == "swine", branchlen <= 15) %>%
    mutate(position = as.numeric(position)) %>%
    group_by(position, residue_to) %>%
    tally() %>%
    dplyr::rename(n_swine_residue = n) %>%
    ungroup()
  
  # Summarize human and swine observed AAs
  human_observed_summary <- left_join(human_observed_residue_tally, total_human_observed_counts, by = c("position")) %>%
    group_by(position) %>%
    arrange(desc(n_human_residue), .by_group = TRUE) %>%
    summarise(
      "Human observed AAs" = paste0(residue_to, "(", round(n_human_residue / total_human_observed * 100, 1), "%)", collapse = ", ")
    ) %>%
    ungroup() %>%
    mutate(segment = segment_name, .before = "position")
  
  swine_observed_summary <- left_join(swine_observed_residue_tally, total_swine_observed_counts, by = c("position")) %>%
    group_by(position) %>%
    arrange(desc(n_swine_residue), .by_group = TRUE) %>%
    summarise(
      "Swine observed AAs" = paste0(residue_to, "(", round(n_swine_residue / total_swine_observed * 100, 1), "%)", collapse = ", ")
    ) %>%
    ungroup() %>%
    mutate(segment = segment_name, .before = "position")
  
  
  # Join the summaries
  branchdiff_summarised <- left_join(residue_summary, mutation_summary, by = c("segment", "position", "trait")) %>%
    left_join(human_observed_summary, by = c("segment","position")) %>%
    left_join(swine_observed_summary, by = c("segment","position"))
  
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

###############

shap_positions_amino_acids <- inner_join(combined_branchdiff_sumarised, positions_all, by = c("segment", "position", "trait"))

write.table(shap_positions_amino_acids, "figures_tables/figure4_5_S4_S5_xgboost/shap_positions_amino_acids_summary.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

