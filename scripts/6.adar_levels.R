library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(matrixStats)


tissues <- list(
  cortex = readRDS('all_sites_cortex_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hippocampus = readRDS('all_sites_hippocampus_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hypothalamus = readRDS('all_sites_hypothalamus_rescue_with_add_filter_g_2_e_005_r_5.rds')
)

overall_editing_level_per_replicate <- map(tissues, function(chrom_list) {
  
  chrom_list %>%
    map(function(rep_list) {
      bind_rows(rep_list, .id = "replicate_id")
    }) %>%
    bind_rows(.id = "chromosome") %>%
    group_by(replicate_id) %>%
    summarise(
      total_G = sum(G_count, na.rm = TRUE),
      total_A = sum(A_count, na.rm = TRUE),
      G_over_A = total_G / (total_A + total_G),
      .groups = "drop"
    )
})

tissue_names <- c('cortex', 'hippocampus', 'hypothalamus')

for (tissue in tissue_names){
  overall_editing_level_per_replicate[[tissue]] <- overall_editing_level_per_replicate[[tissue]] %>%
    mutate(temp_num = as.numeric(gsub("\\D", "", replicate_id))) %>% 
    arrange(temp_num) %>%
    select(-temp_num)
}

vector_overall_editing_levels <- c()

for (tissue in tissue_names){
  vector_overall_editing_levels <- append(vector_overall_editing_levels, 
                                    overall_editing_level_per_replicate[[tissue]]$G_over_A)
}


counts_table <- read.table('Final_counts_CPM_brain.v112.txt')
counts_table <- counts_table %>% relocate(6:8)


adar_ids <- c('ENSSSCG00000006543', 'ENSSSCG00000012089', 'ENSSSCG00000011160')

adar_data_list <- list(
  adar1 = as.numeric(counts_table[adar_ids[1], ]),
  adar2 = as.numeric(counts_table[adar_ids[2], ]),
  adar3 = as.numeric(counts_table[adar_ids[3], ])
)


col_names <- colnames(counts_table)

sorted_adars <- list()

for (i in 1:3) {

  current_adar <- adar_data_list[[i]]
  

  cort <- current_adar[endsWith(col_names, "_C")]
  hipp <- current_adar[endsWith(col_names, "_E")]
  hyp  <- current_adar[!endsWith(col_names, "_C") & !endsWith(col_names, "_E")]
  
  sorted_adars[[paste0("adar", i)]] <- c(cort, hipp, hyp)
}


plot_adar_data <- data.frame(adar1_expression = sorted_adars[[1]],
                        adar2_expression = sorted_adars[[2]],
                        adar3_expression = sorted_adars[[3]],
                        overall_editing_levels= vector_overall_editing_levels)


n_total <- nrow(plot_adar_data)


tissue_vector <- c(
  rep("Cortex", 48),          
  rep("Hippocampus", 48),       
  rep("Hypothalamus", n_total - 96) 
  )

plot_adar_data$Tissue <- factor(tissue_vector, levels = c("Cortex", "Hippocampus", "Hypothalamus"))
number_of_editing_sites_per_rep <- ncol(editing_matrix) - rowCounts(editing_matrix, value = 0) 
plot_adar_data$Number_of_Sites <- number_of_editing_sites_per_rep

first_two <- rep(1:48, times = 2)  
third <- setdiff(1:48, c(2, 10, 14))
third_rep <- third
animal_ids <- c(first_two, third_rep)
plot_adar_data$PigID <- animal_ids
#plot_data$Norm_N_of_sites <- normalized_n_of_sites_per_million



##ADAR levels

clean_data <- tibble(
  ADAR1 = sorted_adars[[1]],
  ADAR2 = sorted_adars[[2]],
  ADAR3 = sorted_adars[[3]],
  Tissue = factor(tissue_vector, levels = c("Cortex", "Hippocampus", "Hypothalamus"))
)


plot_adar <- clean_data %>%
  pivot_longer(
    cols = c(ADAR1, ADAR2, ADAR3),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  group_by(Tissue, Gene) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SD_Expression   = sd(Expression),
    .groups = "drop"
  )

ggplot(plot_adar, aes(x = Tissue, y = Mean_Expression, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = Mean_Expression - SD_Expression,
        ymax = Mean_Expression + SD_Expression),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  labs(
    title = "Mean ADAR Expression by Tissue",
    y = "Mean Expression Levels (TPM)",
    x = "Tissue Type"
  )



##statistics summary

stats_overall_editing_levels <- plot_data_long %>%
  group_by(ADAR, Tissue) %>%
  do(tidy(lm(overall_editing_levels ~ Expression, data = .))) %>%
  filter(term == "Expression") %>%
  select(ADAR, Tissue, estimate, p.value) %>%
  left_join(
    plot_data_long %>%
      group_by(ADAR, Tissue) %>%
      do(glance(lm(overall_editing_levels ~ Expression, data = .))) %>%
      select(ADAR, Tissue, r.squared),
    by = c("ADAR", "Tissue")
  ) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

stats_sites <- plot_data_long %>%
  group_by(ADAR, Tissue) %>%
  do(tidy(lm(Number_of_Sites ~ Expression, data = .))) %>%
  filter(term == "Expression") %>%
  select(ADAR, Tissue, estimate, p.value) %>%
  left_join(
    plot_data_long %>%
      group_by(ADAR, Tissue) %>%
      do(glance(lm(Number_of_Sites ~ Expression, data = .))) %>%
      select(ADAR, Tissue, r.squared),
    by = c("ADAR", "Tissue")
  ) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

print(stats_sites)
print(stats_overall_editing_levels)

write_tsv(stats_overall_editing_levels, 'model_stats_summary_oel.tsv')
