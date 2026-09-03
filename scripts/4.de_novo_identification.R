library(tibble)
library(readr)
library(dplyr)
library(GenomicRanges)
library(JACUSA2helper)
library(tidyr)
library(purrr)

tissues <- c('cortex', 'hippocampus', 'hypothalamus')
chromosomes <- 1:18

# Filters
MIN_COVERAGE_DEPTH <- 10    
MIN_G_COUNT        <- 2     
MIN_RATE           <- 0.05  
MIN_CONSENSUS_REPS <- 5     


master_tracking_list <- list() # this is for capturing the amount of candidate
                               # editing sites that are left after the 
                               # application of each filter step

#Note: the JACUSA output files need to be in the same directory as the script


track_idx <- 1

for (tissue in tissues) {
  
  master_list_sites_with_rates_all_chroms <- list()
  
  for (chrom in chromosomes) {
    
    data_file <- paste0(tissue, '_calls_chrom_', chrom, '.out')
    discarded_file <- paste0(tissue, '_calls_chrom_', chrom, '.out.filtered')
    
   # track step 1: count all the candidate sites in the main jacusa file with
   # an A reference allele
    
    data_raw <- read_result(data_file)
    #cond_name <- if("cond1" %in% names(data_raw$bases)) "cond1" else names(data_raw$bases)[1]
    cond_name <- 'cond1'
    replicate_names <- names(data_raw$bases[[cond_name]])
    
    data <- data_raw[data_raw$ref == "A", ]
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "1_Main_File_Ref_A",
      Sites_Processed = length(data_raw), Sites_Passed = length(data), 
      Sites_Filtered = length(data_raw) - length(data)
    )
    track_idx <- track_idx + 1
    
    # track step 2: count all the candidate sites with more than two alleles
    #present and A reference allele in the discarded jacusa file
    sites_disc_raw <- 0
    sites_disc_filt <- 0
    
    
    discarded_raw <- read_result(discarded_file)
    sites_disc_raw <- length(discarded_raw)
      
    discarded <- discarded_raw[discarded_raw$filter == 'M' & discarded_raw$ref == 'A']
    sites_disc_filt <- length(discarded)
      
    combined_data <- c(data, discarded)
    combined_data <- data
    
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "2_Discarded_File_M_A",
      Sites_Processed = sites_disc_raw, Sites_Passed = sites_disc_filt, 
      Sites_Filtered = sites_disc_raw - sites_disc_filt
    )
    track_idx <- track_idx + 1
    
    combined_data <- sort(combined_data)
    sites_combined <- length(combined_data)
    
    cov_matrix <- sapply(combined_data$bases[[cond_name]], rowSums)
    
    #to shorten the computational time, keep only sites where at least one
    #replicate had coverage greater or equal to 10. Wasn't really that useful
    
    keep_indices <- rowSums(cov_matrix >= MIN_COVERAGE_DEPTH) >= 1
    cov_filtered_data <- combined_data[keep_indices, ]
    
    sites_cov_filt <- length(cov_filtered_data)
    
    
    #track step 3: count candidates with at least one replicate having coverage
    #more than 10
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "3_Min_Coverage_Any_Rep",
      Sites_Processed = sites_combined, Sites_Passed = sites_cov_filt, 
      Sites_Filtered = sites_combined - sites_cov_filt
    )
    track_idx <- track_idx + 1

    
    # phase 1
    print(paste("Chrom", chrom, ": Phase 1 Discovery"))
    
    a_to_g_ids_list <- list()
    candidate_ids_list <- list()
    
    for (rep in replicate_names) {
      
      current_rep_bases <- cov_filtered_data$bases[[cond_name]][[rep]]
      rep_coverage <- rowSums(as.matrix(current_rep_bases))
      cov_pass <- rep_coverage >= MIN_COVERAGE_DEPTH
      
      if(sum(cov_pass) == 0) next
      
      base_calls <- base_call(current_rep_bases)
      base_subs  <- base_sub(base_calls, cov_filtered_data$ref)
      is_a_to_g  <- !is.na(base_subs) & (base_subs == "A->G")
      
      discovery_indices <- cov_pass & is_a_to_g
      
      if(sum(discovery_indices) > 0) {
        a_to_g_ids_list[[rep]] <- JACUSA2helper::id(cov_filtered_data[discovery_indices])
      }
      
      if(sum(discovery_indices) == 0) next
      
      temp_sites <- cov_filtered_data[discovery_indices]
      temp_matrix <- as.matrix(temp_sites$bases[[cond_name]][[rep]])
      
      G_counts <- temp_matrix[, "G"]
      Total_counts <- rowSums(temp_matrix)
      Rates <- G_counts / Total_counts
      
      strict_pass <- (G_counts >= MIN_G_COUNT) & (Rates >= MIN_RATE)
      
      if(sum(strict_pass) > 0) {
        candidate_ids_list[[rep]] <- JACUSA2helper::id(temp_sites[strict_pass])
      }
    }
    
  
    # track step 4: count candidate sites with A-to-G(I) substitution
    
    all_a_to_g_ids <- unlist(a_to_g_ids_list)
    sites_a_to_g <- length(unique(all_a_to_g_ids))
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "4_Phase1_A_to_G_Only",
      Sites_Processed = sites_cov_filt, Sites_Passed = sites_a_to_g, 
      Sites_Filtered = sites_cov_filt - sites_a_to_g
    )
    track_idx <- track_idx + 1
    
    
    #track step 5: count candidate sites with editing rate higher than 0.05
    #and at least two edited reads
    
    all_passing_ids <- unlist(candidate_ids_list)
    sites_strict <- length(unique(all_passing_ids))
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "5_Phase1_Strict_Thresholds",
      Sites_Processed = sites_a_to_g, Sites_Passed = sites_strict, 
      Sites_Filtered = sites_a_to_g - sites_strict
    )
    track_idx <- track_idx + 1
    
    #phase 2
    
    print("Phase 2: Consensus Check")
    
    id_counts <- table(all_passing_ids)
    valid_site_ids <- names(id_counts[id_counts >= MIN_CONSENSUS_REPS])
    
    # tracking step 6: count sites where at least 5 tissue replicates satisfy
    # the minimum editing rate and minimum number of edited reads criteria
    
    sites_consensus <- length(valid_site_ids)
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "6_Phase2_Consensus",
      Sites_Processed = sites_strict, Sites_Passed = sites_consensus, 
      Sites_Filtered = sites_strict - sites_consensus
    )
    track_idx <- track_idx + 1
    
    print(paste("Initial Valid Sites:", length(valid_site_ids)))
    if(length(valid_site_ids) == 0) next
    
    
    #phase 3
    print("Phase 3: Data Retrieval")
    
    final_sites_object <- cov_filtered_data[JACUSA2helper::id(cov_filtered_data) %in% valid_site_ids]
    
    sites_before_snp <- length(final_sites_object)
    
    n_sites <- length(final_sites_object)
    n_reps  <- length(replicate_names)
    
    rate_matrix <- matrix(NA, nrow = n_sites, ncol = n_reps)
    colnames(rate_matrix) <- replicate_names
    
    for (i in seq_along(replicate_names)) {
      rep_name <- replicate_names[i]
      mat <- as.matrix(final_sites_object$bases[[cond_name]][[rep_name]])
      
      g <- mat[, "G"]
      tot <- rowSums(mat)
      
      r <- g / tot
      r[tot == 0] <- NA 
      rate_matrix[, i] <- r
    }
    
    has_full_editing <- rowSums(rate_matrix > 0.999, na.rm = TRUE) > 0
    
    rate_matrix_nonzero <- rate_matrix
    rate_matrix_nonzero[rate_matrix_nonzero == 0] <- NA
    mean_rates_nonzero <- rowMeans(rate_matrix_nonzero, na.rm = TRUE)
    is_high_mean <- mean_rates_nonzero > 0.95
    is_high_mean[is.na(is_high_mean)] <- FALSE
    
    snp_suspects <- has_full_editing & is_high_mean
    
    if (sum(snp_suspects) > 0) {
      print(paste("Removing", sum(snp_suspects), "sites likely to be genomic variants."))
      final_sites_object <- final_sites_object[!snp_suspects]
      valid_site_ids <- JACUSA2helper::id(final_sites_object)
    }
    
    #track 7: count detected editing sites after removal of sites likely
    #to be genomic variants 
    sites_after_snp <- length(final_sites_object)
    
    master_tracking_list[[track_idx]] <- tibble(
      Tissue = tissue, Chromosome = chrom, Step = "7_Phase3_SNP_Filter",
      Sites_Processed = sites_before_snp, Sites_Passed = sites_after_snp, 
      Sites_Filtered = sites_before_snp - sites_after_snp
    )
    track_idx <- track_idx + 1
    
    replicate_sites_and_rates_for_chrom <- list()
    
    for (rep in replicate_names) {
      
      rep_matrix <- as.matrix(final_sites_object$bases[[cond_name]][[rep]])
      
      rep_coverage <- rowSums(rep_matrix)
      rep_g_counts <- rep_matrix[, "G"]
      
      rescue_pass <- (rep_coverage >= MIN_COVERAGE_DEPTH) & (rep_g_counts >= 1)
      
      if(sum(rescue_pass) == 0) {
        empty_tibble <- tibble(
          site_id = character(), 
          editing_rate = numeric(), 
          A_count = numeric(), 
          G_count = numeric(), 
          coverage = numeric()
        )
        replicate_sites_and_rates_for_chrom[[rep]] <- empty_tibble
        next 
      }
      
      clean_matrix <- rep_matrix[rescue_pass, , drop=FALSE]
      
      A_counts <- clean_matrix[, "A"]
      G_counts <- clean_matrix[, "G"]
      Total_counts <- rowSums(clean_matrix)
      Rates <- G_counts / Total_counts
      
      subset_object <- final_sites_object[rescue_pass]
      
      summary_tibble <- tibble(
        site_id      = JACUSA2helper::id(subset_object),
        editing_rate = Rates,
        A_count      = A_counts,
        G_count      = G_counts,
        coverage     = Total_counts
      )
      
      replicate_sites_and_rates_for_chrom[[rep]] <- summary_tibble
    }
    
    
    list_name <- paste0("chrom_", chrom)
    master_list_sites_with_rates_all_chroms[[list_name]] <- replicate_sites_and_rates_for_chrom
    
    text_file_name <- paste0(tissue, '_chrom_', chrom, '_test_for_filter_counts_two.txt')
    readr::write_lines(valid_site_ids, text_file_name)
    
    print("--- Chromosome processed ---")
  }
  
  output_filename <- paste0('all_sites_', tissue, '_test_for_filter_counts.rds')
  saveRDS(master_list_sites_with_rates_all_chroms, file = output_filename)
}


#summarize the filtering steps

filter_summary_df <- bind_rows(master_tracking_list)
write_csv(filter_summary_df, "RNA_editing_filter_summary.csv")
print("Filtering summary saved to RNA_editing_filter_summary.csv")

summary <- read_csv("RNA_editing_filter_summary.csv")

aggregated_summary <- summary %>%
  group_by(Tissue, Step) %>%
  summarise(
    Total_Sites_Surviving = sum(Sites_Passed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Tissue, Step)

print(aggregated_summary)

sites_wide <- aggregated_summary %>%
  select(Tissue, Step, Total_Sites_Surviving) %>%
  pivot_wider(names_from = Tissue, values_from = Total_Sites_Surviving)

cat("\n--- Surviving Sites per Tissue ---\n")
print(sites_wide)

#for retrieving the Q607R editing sites in GRIA2 to be included later on

all_sites_snps <- snps %>%
    flatten() %>% 
    keep(~ is.data.frame(.x) && "site_id" %in% names(.x)) %>%
    map("site_id") %>%
    unlist(recursive = TRUE)
  
unique_sites_snps <- unique(all_sites_snps)
writeLines(unique_sites_snps, "unique__site_snps_ids.txt")
