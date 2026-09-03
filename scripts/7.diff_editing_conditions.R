library(broom)
library(dplyr)
library(purrr)
library(readxl)
library(stringr)

meta_raw <- read_excel("Info_samples_brain.xlsx")

meta <- meta_raw %>%
  mutate(
    Pig_ID = as.integer(Pig_ID),
    Tissue = case_when(
      str_detect(Tissue, regex("^hypothalamus", ignore_case = TRUE)) ~ "hypothalamus",
      str_detect(Tissue, regex("^cortex", ignore_case = TRUE)) ~ "cortex",
      str_detect(Tissue, regex("^hippocampus", ignore_case = TRUE)) ~ "hippocampus",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    Sample_Name = `Sample Name`,
    Pig_ID,
    Tissue,
    Treatment,
    Allostatic_load
  )

table(meta$Tissue)

matrix_meta <- tibble(
  sample_id = rownames(editing_matrix),
  Tissue = str_extract(sample_id, "(cortex|hippocampus|hypothalamus)")
)

pig_map <- meta %>%
  arrange(Tissue, Pig_ID) %>%
  group_by(Tissue) %>%
  mutate(tissue_index = row_number()) %>%
  ungroup()

matrix_meta <- matrix_meta %>%
  group_by(Tissue) %>%
  mutate(tissue_index = row_number()) %>%
  ungroup()

final_metadata <- matrix_meta %>%
  left_join(pig_map, by = c("Tissue", "tissue_index")) %>%
  select(
    sample_id,
    Pig_ID,
    Tissue,
    Treatment,
    Allostatic_load
  )

run_condition_mw <- function(tissue, condition_col, edit_mat, metadata) {
  
  tissue_samples <- metadata %>%
    filter(.data$Tissue == tissue)
  
  conditions <- sort(unique(tissue_samples[[condition_col]]))
  
  s1 <- tissue_samples %>% filter(.data[[condition_col]] == conditions[1])
  s2 <- tissue_samples %>% filter(.data[[condition_col]] == conditions[2])
  
  
  mat1 <- edit_mat[s1$sample_id, , drop = FALSE]
  mat2 <- edit_mat[s2$sample_id, , drop = FALSE]
  
  
  res <- purrr::map_dfr(colnames(edit_mat), function(site) {
    
    x <- as.numeric(mat1[, site])
    y <- as.numeric(mat2[, site])
    
    
    # if a site is 0 in almost all editing events it's likely noise
    n_nonzero <- sum(c(x, y) > 0)
    
    
    # at least 3 edited samples
    if (n_nonzero < 3) {
      return(tibble(
        site = site,
        p_value = NA_real_,
        n_group1 = length(x),
        n_group2 = length(y),
        median_diff = NA_real_
      ))
    }
    
    
    if (length(unique(c(x, y))) < 2) {
      return(tibble(
        site = site,
        p_value = NA_real_,
        n_group1 = length(x),
        n_group2 = length(y),
        median_diff = NA_real_
      ))
    }
    
    
    wt <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
    
    
    tibble(
      site = site,
      p_value = wt$p.value,
      n_group1 = length(x),
      n_group2 = length(y),
      median_diff = median(x) - median(y)
    )
  })
  
  
  res %>%
    mutate(
      tissue = tissue,
      condition_comparison = paste(conditions[1], "vs", conditions[2]),
      p_adj = p.adjust(p_value, method = "BH")
    )
}



test_tissues <- c("cortex", "hippocampus", "hypothalamus")

mann_whitney_results <- purrr::map_dfr(
  test_tissues,
  ~run_condition_mw(.x, "Treatment", editing_matrix, final_metadata)
)

print(mann_whitney_results %>%
        filter(p_adj< 0.05)) 


run_glm_per_tissue <- function(tissue_name, edit_mat, metadata,
                               min_nonzero = 3) {
  
  tissue_meta <- metadata %>% filter(Tissue == tissue_name)
  
  tissue_meta$Treatment <- factor(tissue_meta$Treatment)
  tissue_meta$Allostatic_load <- factor(tissue_meta$Allostatic_load)
  
  sample_ids <- tissue_meta$sample_id
  
  map_dfr(colnames(edit_mat), function(site) {
    
    y <- as.numeric(edit_mat[sample_ids, site])
    keep <- complete.cases(y)
    
    y <- y[keep]
    meta_sub <- tissue_meta[keep, ]
    

    fit <- glm(y ~ Treatment * Allostatic_load, data = meta_sub)
    
    tf <- broom::tidy(fit)
    
    tibble(
      site = site,
      treatment_p = tf$p.value[grepl("^Treatment", tf$term)],
      allostatic_p = tf$p.value[grepl("^Allostatic_load", tf$term)],
      interaction_p = tf$p.value[grepl(":", tf$term)]
    )
  }) %>%
    mutate(
      tissue = tissue_name,
      treatment_p_adj = p.adjust(treatment_p, "BH"),
      allostatic_p_adj = p.adjust(allostatic_p, "BH"),
      interaction_p_adj = p.adjust(interaction_p, "BH")
    )
}

tissue_names <- unique(final_metadata$Tissue)

all_glm_results <- purrr::map_dfr(
  tissue_names,
  ~run_glm_per_tissue(.x, editing_matrix, final_metadata)
)

all_glm_results %>%
  filter(allostatic_p_adj < 0.05)

run_glm_per_tissue_single_factor <- function(tissue_name,
                                             edit_mat,
                                             metadata,
                                             factor_var,
                                             min_nonzero = 3) {
  
  tissue_meta <- metadata %>%
    filter(Tissue == tissue_name)
  
  tissue_meta[[factor_var]] <- factor(tissue_meta[[factor_var]])
  
  sample_ids <- tissue_meta$sample_id
  
  map_dfr(colnames(edit_mat), function(site) {
    
    y <- as.numeric(edit_mat[sample_ids, site])
    keep <- complete.cases(y)
    
    y <- y[keep]
    meta_sub <- tissue_meta[keep, ]
    
    
    form <- as.formula(paste("y ~", factor_var))
    
    
    fit <- glm(form, data = meta_sub)
    
    tf <- broom::tidy(fit)
    
    tibble(
      site = site,
      factor = factor_var,
      p_value = tf$p.value[tf$term != "(Intercept)"]
    )
  }) %>%
    mutate(
      tissue = tissue_name,
      p_adj = p.adjust(p_value, method = "BH")
    )
}

tissue_names <- unique(final_metadata$Tissue)

all_glm_results_single <- purrr::map_dfr(
  tissue_names,
  ~run_glm_per_tissue_single_factor(
    tissue_name = .x,
    edit_mat = editing_matrix,
    metadata = final_metadata,
    factor_var = "Treatment"
  )
)

all_glm_results_single %>%
  filter(p_adj < 0.05)