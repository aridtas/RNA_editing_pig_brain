library(tidyverse)


tissues <- list(
  cortex = readRDS('all_sites_cortex_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hippocampus = readRDS('all_sites_hippocampus_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hypothalamus = readRDS('all_sites_hypothalamus_rescue_with_add_filter_g_2_e_005_r_5.rds')
)

##matrix construction


long_df <- tissues %>%
  map(function(chrom_list) {
    
    chrom_list %>%
      map(function(rep_list) {
        
        bind_rows(rep_list, .id = "replicate_id")
        
      }) %>%
      
      list_rbind(names_to = "chromosome")
    
  }) %>%
  
  list_rbind(names_to = "tissue")

long_df <- long_df %>%
  mutate(unique_sample_id = paste(tissue, replicate_id, sep = "_"))

print(paste("Total rows in the matrix:", nrow(long_df)))

wide_df <- long_df %>%
  dplyr::select(unique_sample_id, site_id, editing_rate) %>%
  pivot_wider(
    names_from = site_id, 
    values_from = editing_rate,
    values_fill = 0 
  )

editing_matrix_data <- wide_df %>% 
  column_to_rownames("unique_sample_id")

editing_matrix <- as.matrix(editing_matrix_data)

non_zero_count <- sum(editing_matrix_data != 0)

count_plus <- sum(grepl("\\+$", colnames(editing_matrix_data)))
count_minus <- sum(grepl("-$", colnames(editing_matrix_data)))


# unique sites

df_unique_sites <- data.frame(colnames(editing_matrix_data))

df_unique_sites_split <- df_unique_sites %>%
  separate(colnames.editing_matrix_data., 
           into = c("Chromosome", "Coordinate"), 
           sep = ":|-",           
           extra = "drop")        

write.table(df_unique_sites_split,'unique_editing_sites.tsv',
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')


##highly edited sites


cortex_idx <- 1:48
hippo_idx <- 49:96
hypo_idx <- 97:141

#function to calculate mean (excluding zeros) and sample count
calculate_mean_and_count <- function(sub_matrix) {
  # Calculate mean of non-zero entries per column
  means <- apply(sub_matrix, 2, function(x) {
    vals <- x[x > 0]
    if(length(vals) == 0) return(0) else return(mean(vals))
  })
  
  # Calculate number of non-zero samples per column
  counts <- colSums(sub_matrix > 0)
  
  list(mean = means, count = counts)
}


cortex_stats <- calculate_mean_and_count(editing_matrix_data[cortex_idx, ])
hippo_stats  <- calculate_mean_and_count(editing_matrix_data[hippo_idx, ])
hypo_stats   <- calculate_mean_and_count(editing_matrix_data[hypo_idx, ])


means_of_editing_sites_per_tissue <- data.frame(
  site = colnames(editing_matrix_data),
  mean_cortex = cortex_stats$mean,
  mean_hippo  = hippo_stats$mean,
  mean_hypo   = hypo_stats$mean,
  count_cortex = cortex_stats$count,
  count_hippo  = hippo_stats$count,
  count_hypo   = hypo_stats$count,
  stringsAsFactors = FALSE
)

highly_edited_sites <- means_of_editing_sites_per_tissue[means_of_editing_sites_per_tissue$mean_cortex > 0.6 | 
                      means_of_editing_sites_per_tissue$mean_hippo  > 0.6 | 
                        means_of_editing_sites_per_tissue$mean_hypo   > 0.6, ]

# reset row names to site names for cleanliness
rownames(highly_edited_sites) <- highly_edited_sites$site
highly_edited_sites$site <- NULL

write.table(rownames(highly_edited_sites),'highly_edited_sites_06.tsv',
            row.names = FALSE,
            quote = FALSE,
            col.names = FALSE,
            sep = '\t')



#the paper constructing a database of pig RNA editing sites
presdb <- read.table('presdb.tsv', sep = '\t', header = TRUE)
presdb <- presdb %>% unite("id", 'Chromsome', 'Position', sep = "_")
unique_edit <- read.table('unique_editing_sites.tsv', sep = '\t', header = TRUE)
unique_edit <- unique_edit %>% unite('id', 'Chromosome', 'Coordinate', sep = "_")
common_elements <- intersect(presdb$id, unique_edit$id)

#the Nature paper
atoi <- read_tsv('A-to-I.tsv')
atoi <- atoi %>% unite('id', 'Chromosome', 'Coordinate', sep = '_')
common_elements2 <- intersect(presdb$id, atoi$id)
common_elements3 <- intersect(unique_edit$id, atoi$id)


