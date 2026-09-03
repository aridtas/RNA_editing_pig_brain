library(tidyverse)
library(ggplot2)

metadata <- tibble(sample_id = rownames(editing_matrix)) %>%
  separate(sample_id,
           into = c("tissue", "pig_id"),
           sep = "_",
           remove = FALSE) %>%
  mutate(
    tissue = factor(tissue),
    pig_id = factor(pig_id)
  )


tissue_names <- c('cortex', 'hippocampus', 'hypothalamus')



sites_counts_per_tissue <- data.frame()

for (tis in tissue_names) {
  

  target_samples <- metadata %>% 
    filter(tissue == tis) %>% 
    pull(sample_id)
  
  
  # subset the matrix
  # intersect is to make sure samples are grabbed that actually exist in the matrix
  sub_matrix <- editing_matrix[target_samples, , drop = FALSE]
  

  site_support_counts <- colSums(sub_matrix > 0, na.rm = TRUE)

  n_samples_in_tissue <- length(target_samples)
  
  # Specific Thresholds
  stats <- list(
    Tissue = tis,
    Total_Samples = n_samples_in_tissue,
    Count_Absent = sum(site_support_counts == 0),
    Count_Edited = sum(site_support_counts > 0),
    # Support Thresholds - wanted to check how many sites are left by changing
    #minimum number of samples
    Support_GE_6  = sum(site_support_counts >= 6),
    Support_GE_12 = sum(site_support_counts >= 12),
    Support_GE_20 = sum(site_support_counts >= 20),
    Support_GE_30 = sum(site_support_counts >= 30),
    Support_GE_40 = sum(site_support_counts >= 40),
    Present_In_ALL = sum(site_support_counts == n_samples_in_tissue)
  )

  sites_counts_per_tissue <- rbind(sites_counts_per_tissue,
                                       as.data.frame(stats))
}


##differential editing between tissues

tissue_pairs <- list(
  c("cortex", "hippocampus"),
  c("cortex", "hypothalamus"),
  c("hippocampus", "hypothalamus")
)

run_wilcoxon_signed_rank <- function(tissue1, tissue2, edit_mat, metadata) {
  

  s1 <- metadata %>% filter(tissue == tissue1)
  s2 <- metadata %>% filter(tissue == tissue2)
  
  #these steps are to account for the three missing hypothalamus samples
  paired_ids <- intersect(s1$pig_id, s2$pig_id)
  s1 <- s1 %>% filter(pig_id %in% paired_ids) #%>% arrange(pig_id)
  s2 <- s2 %>% filter(pig_id %in% paired_ids) #%>% arrange(pig_id)
  

  mat1 <- edit_mat[s1$sample_id, , drop = FALSE]
  mat2 <- edit_mat[s2$sample_id, , drop = FALSE]
  

  
  res <- purrr::map_dfr(colnames(edit_mat), function(site) {
    
    # extract vectors for this site
    x <- as.numeric(mat1[, site]) 
    y <- as.numeric(mat2[, site]) 
    
    # drop NAs
    keep <- complete.cases(x, y)
    x <- x[keep]
    y <- y[keep]
    
    #at least 12 samples per tissue for site to be ok for testing
    
    count_x_nonzero <- sum(x > 0, na.rm = TRUE)
    count_y_nonzero <- sum(y > 0, na.rm = TRUE)
    
    if (count_x_nonzero < 12 | count_y_nonzero < 12) {
      # if not enough data, return NA row and skip the testing
      return(tibble(
        site = site,
        p_value = NA_real_,
        n_pairs = length(x),
        median_diff = NA_real_,
        status = "Filtered_Low_Count"
      ))
    }
    
    
    wt <- wilcox.test(x, y, paired = TRUE, exact = FALSE)
  
    
    tibble(
      site = site,
      p_value = wt$p.value,
      n_pairs = length(x),
      median_diff = median(x - y), # Magnitude of change
      status = "Tested"
    )
  })
  
  # BH multiple testing correction
  
  res <- res %>%
    mutate(
      tissue_comparison = paste(tissue1, "vs", tissue2),
      p_adj = p.adjust(p_value, method = "BH") 
    )
  
  return(res)
}


results_wilcoxon <- map_dfr(
  tissue_pairs,
  ~run_wilcoxon_signed_rank(.x[1], .x[2], editing_matrix, metadata)
)

significant_changes <- results_wilcoxon %>% 
  filter(p_adj < 0.05) %>%
  arrange(p_adj)


significant_counts <- significant_changes %>%
  dplyr::count(tissue_comparison)
print(significant_counts)


sig_cort_hipp <- significant_changes %>% 
  filter(tissue_comparison == "cortex vs hippocampus") %>% 
  pull(site)

sig_cort_hypo <- significant_changes %>% 
  filter(tissue_comparison == "cortex vs hypothalamus") %>% 
  pull(site)

sig_hipp_hypo <- significant_changes %>% 
  filter(tissue_comparison == "hippocampus vs hypothalamus") %>% 
  pull(site)


shared_all_three <- Reduce(intersect, list(sig_cort_hipp, sig_cort_hypo, 
                                           sig_hipp_hypo))
shared_cort_hypo_AND_hipp_hypo <- intersect(sig_cort_hypo, sig_hipp_hypo)



results_wilcoxon %>%
  group_by(tissue_comparison) %>%
  summarise(
    tested_sites = sum(!is.na(p_value)),
    total_sites = n()
  )

results_wilcoxon %>%
  filter(!is.na(p_adj)) %>%
  group_by(tissue_comparison) %>%
  summarise(
    frac_significant = mean(p_adj < 0.05)
  )

cort_v_hipp <- results_wilcoxon %>%
  filter(p_adj< 0.05) %>%
  filter(tissue_comparison == "cortex vs hippocampus")

cort_v_hyp <- results_wilcoxon %>%
  filter(p_adj< 0.05) %>%
  filter(tissue_comparison == "cortex vs hypothalamus")


hipp_v_hyp <- results_wilcoxon %>%
  filter(p_adj< 0.05) %>%
  filter(tissue_comparison == "hippocampus vs hypothalamus")

readr::write_lines(cort_v_hipp[[1]], "cort_v_hipp_wilcoxon.txt")
readr::write_lines(cort_v_hyp[[1]], "cort_v_hyp_wilcoxon.txt")
readr::write_lines(hipp_v_hyp[[1]], "hipp_v_hyp.txt")


#Volcano Plot


# Create a plotting data frame from the FULL results
# We add a 'sig' column to label points as Significant/Not Significant
volcano_data <- results_wilcoxon %>%
  # Remove NAs (sites that were filtered out by the 12-sample rule)
  filter(!is.na(p_adj)) %>%
  mutate(
    sig = ifelse(p_adj < 0.05, "Significant", "Not Significant"),)



ggplot(volcano_data, aes(x = median_diff, y = -log10(p_adj))) +
  geom_point(data = subset(volcano_data, sig == "Not Significant"), 
             color = "grey70", alpha = 0.5, size = 1) +
  geom_point(data = subset(volcano_data, sig == "Significant"), 
             color = "red", alpha = 0.8, size = 1.5) +
  
  facet_wrap(~ tissue_comparison) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  labs(
    title = "Differential Editing Between Tissues",
    x = "Δ Editing Rate (Median Difference)",
    y = "-log10(BH-adjusted p-value)"
  ) +
  theme_bw()

write.table(unique(significant_changes$site),
            'diff_edited_sites.tsv',
            row.names = FALSE,
            col.names = FALSE,
            sep ='\t',
            quote = FALSE)

cort_v_hyp <- significant_changes%>%
  filter(tissue_comparison == "cortex vs hypothalamus")

readr::write_lines(cort_v_hyp[[1]], "cort_v_hipp_wilcoxon.txt")


tissue_pairs <- list(
  c("cortex", "hippocampus"),
  c("cortex", "hypothalamus"),
  c("hippocampus", "hypothalamus")
)


# Count the number of positive and negative median differences per comparison
direction_counts <- significant_changes %>%
  mutate(direction = ifelse(median_diff > 0, "Tissue 1 > Tissue 2",
                            "Tissue 1 < Tissue 2")) %>%
  group_by(tissue_comparison, direction) %>%
  summarise(count = n(), .groups = "drop")

print(direction_counts)

