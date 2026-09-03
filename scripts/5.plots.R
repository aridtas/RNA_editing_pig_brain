library(purrr)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(readr)
library(rstatix)
library(tidyr)
library(matrixStats)
library(tibble)
library(scales)

tissues <- list(
  cortex = readRDS('all_sites_cortex_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hippocampus = readRDS('all_sites_hippocampus_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hypothalamus = readRDS('all_sites_hypothalamus_rescue_with_add_filter_g_2_e_005_r_5.rds')
)

tissue_names <- c('cortex', 'hippocampus', 'hypothalamus')

##venn diagram

editing_sites_per_tissue <- list()

for (tissue in names(tissues)) {
  
  unique_sites <- tissues[[tissue]] %>%
    map(bind_rows) %>%   
    bind_rows() %>%      
    distinct(site_id) %>%
    pull(site_id)
  
  editing_sites_per_tissue[[tissue]] <- unique_sites
}

venn.diagram(editing_sites_per_tissue, filename = 'venn_diagram_with_add_filter_g_2_e_005_r_5.tiff',
              title ='Unique and shared editing sites across tissues' )

##editing metrics

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

for (tissue in names(tissues)){
  overall_editing_level_per_replicate[[tissue]] <- overall_editing_level_per_replicate[[tissue]] %>%
    mutate(temp_num = as.numeric(gsub("\\D", "", replicate_id))) %>% 
    arrange(temp_num) %>%
    dplyr::select(-temp_num)
}


editing_metrics_plot_df <- bind_rows(
  overall_editing_level_per_replicate,
  .id = "tissue"
)



##normalized number of edited sites by unique mapped reads


df_bam_unique_read_counts <- read.table('bam_read_counts.tsv', sep= '\t', header = TRUE)

cort_read_counts <- df_bam_unique_read_counts[endsWith(df_bam_unique_read_counts$Sample_Name, "_C"), ]
hipp_read_counts <- df_bam_unique_read_counts[endsWith(df_bam_unique_read_counts$Sample_Name, "_E"), ]
hyp_read_counts <- df_bam_unique_read_counts %>% 
  filter(!endsWith(Sample_Name, "_C"), 
         !endsWith(Sample_Name, "_E"))

ordered_read_counts <- rbind(cort_read_counts, hipp_read_counts, hyp_read_counts)

number_of_editing_sites_per_rep <- ncol(editing_matrix) - rowCounts(editing_matrix, value = 0)
normalized_n_of_sites <- number_of_editing_sites_per_rep / ordered_read_counts$Read_Count
normalized_n_of_sites_per_million <- number_of_editing_sites_per_rep / ordered_read_counts$Read_Count * 1000000




number_of_editing_sites_per_rep <- ncol(editing_matrix) - rowCounts(editing_matrix, value = 0)
editing_metrics_plot_df$n_of_sites <- number_of_editing_sites_per_rep




ggplot(editing_metrics_plot_df, aes(x = tissue, y = G_over_A)) +
  geom_violin(trim = FALSE, fill = "grey") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic() +
  labs(
    x = "Distributions of samples per tissue",
    y = "Overall editing level",
    title = "Per-sample Overall editing levels across tissues"
  )

ggplot(editing_metrics_plot_df, aes(x = tissue, y = n_of_sites)) +
  geom_violin(trim = FALSE, fill = "grey") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic() +
  labs(
    x = "Distributions of samples per tissue",
    y = "Number of editing sites",
    title = "Per-sample number of editing sites across tissues"
  )


## friedman

editing_metrics_plot_df$norm_n_of_sites <- normalized_n_of_sites_per_million


# Identify which replicate_ids (pigs) are present in ALL 3 tissues
valid_replicates <- editing_metrics_plot_df %>%
  group_by(replicate_id) %>%
  tally() %>%            # count how many rows each ID has
  filter(n == 3) %>%     # keep only IDs that appear exactly 3 times
  pull(replicate_id)

# filter dataframe to keep only the pigs that have samples in all three tissues,
#this is necessary for friedman test to work

balanced_data <- editing_metrics_plot_df %>%
  filter(replicate_id %in% valid_replicates) %>%
  arrange(replicate_id, tissue) 


run_global_stats <- function(data, measure_var) {
  

  res_friedman <- data %>%
    friedman_test(as.formula(paste(measure_var, "~ tissue | replicate_id")))
  
  print(res_friedman)
  
  # if friedman is significant, do pairwise wilcoxon next
  
  if (res_friedman$p < 0.05) {
    
    pwc <- data %>%
      pairwise_wilcox_test(
        as.formula(paste(measure_var, "~ tissue")),
        paired = TRUE,
        p.adjust.method = "BH"
      )
    
    print(pwc)
    return(pwc)
    
  } else {
    print("friedman wasnt significant")
    return(NULL)
  }
}


res_level  <- run_global_stats(balanced_data, "G_over_A")
res_counts <- run_global_stats(balanced_data, "n_of_sites")
res_normal_counts <- run_global_stats(balanced_data, "norm_n_of_sites" )


#did not include these metrics in the report
balanced_data %>%
  friedman_test(G_over_A ~ tissue | replicate_id) %>%
  mutate(
    kendalls_W = statistic / 
      (n() * (length(unique(balanced_data$tissue)) - 1))
  )

wilcox_effsize(
  balanced_data,
  G_over_A ~ tissue,
  paired = TRUE
)


wilcox_effsize(
  balanced_data,
  n_of_sites ~ tissue,
  paired = TRUE
)


wilcox_effsize(
  balanced_data,
  norm_n_of_sites ~ tissue,
  paired = TRUE
)


## pie chart syn and non-syn mutations

df_mutations <- data.frame(
  category = c("Synonymous", "Non-synonymous"),
  count = c(125, 207)
)

df_mutations <- df_mutations %>%
  mutate(
    percent = count / sum(count),
    label = paste0(
      category, "\n",
      count, " (", scales::percent(percent, accuracy = 1), ")"
    )
  )

ggplot(df_mutations, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  scale_fill_manual(
    values = c("Synonymous" = "#d95f02", "Non-synonymous" = "#1b9e77")
  ) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "Category",
       title = 'Consequences of editing sites found in CDS regions')  +
  theme(plot.title = element_text(size = 15))
  



## aa_changes

aa_changes_files <- c("aa_changes_cortex.tsv",
                      "aa_changes_hypothalamus.tsv",
                      "aa_changes_hippocampus.tsv")

aa_changes_df_list <- lapply(aa_changes_files, function(f) {
  read.table(f, sep = "\t", stringsAsFactors = FALSE)
})

aa_changes_df_all <- do.call(rbind, aa_changes_df_list)


# deduplicate based on V1
aa_changes_unique <- aa_changes_df_all %>%
  distinct(V1, .keep_all = TRUE)

print(paste("Unique sites:", nrow(aa_changes_unique)))


aa_changes_count_data <- aa_changes_unique %>%
  dplyr::count(V2, name = "Count") %>%  
  arrange(desc(Count))                  


colnames(aa_changes_count_data)[1] <- "AA_Change"

ggplot(
  aa_changes_count_data,
  aes(x = reorder(AA_Change, -Count), y = Count)
) +
  geom_col(fill = "#1b9e77") +
  labs(
    title = "Amino Acid changes of unique recoding sites across tissues",
    x = "Amino Acid Change",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


##genomic regions


file_list <- list.files(pattern = "\\.bed$")
file_list <- file_list[!grepl("repetitive", file_list)]

hierarchy_ranks <- c(
  "cds"        = 1,
  "utr3"       = 2, 
  "utr5"       = 3, 
  "introns"    = 4, 
  "intergenic" = 5
)

read_and_label <- function(filename) {
  
  tissue_name <- str_split(filename, "\\.")[[1]][1]
  
  region_name <- dplyr::case_when(
    grepl("cds", filename, ignore.case = TRUE) ~ "cds",
    grepl("utr3", filename, ignore.case = TRUE) ~ "utr3",
    grepl("utr5", filename, ignore.case = TRUE) ~ "utr5",
    grepl("intron", filename, ignore.case = TRUE) ~ "introns",
    grepl("intergenic", filename, ignore.case = TRUE) ~ "intergenic",
    TRUE ~ "unknown"
  )
  
  if (region_name == "unknown") {
    warning(paste("Skipping file with unknown region:", filename))
    return(tibble()) 
  }
  
  df <- read_tsv(filename, col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::select(X1, X2, X3) %>%
    dplyr::rename(chrom = X1, start = X2, end = X3) %>%
    mutate(
      tissue = tissue_name,
      region = region_name,
      rank   = hierarchy_ranks[region_name]
    )
  
  return(df)
}


all_sites_raw <- map_dfr(file_list, read_and_label)

unique_sites_annotated <- all_sites_raw %>%
  mutate(site_id = paste(chrom, start, sep = "_")) %>%
  group_by(site_id) %>%
  slice_min(order_by = rank, n = 1, with_ties = FALSE) %>%
  ungroup()


print(table(unique_sites_annotated$region))

gen_regions_plot_data <- unique_sites_annotated %>%
  dplyr::count(region) %>%
  mutate(region = factor(region, levels = c("cds", "utr5", "utr3", "introns", "intergenic"))) %>%
  mutate(percentage = n / sum(n) * 100)

p <- ggplot(gen_regions_plot_data, aes(x = region, y = n)) +
  
  geom_bar(stat = "identity", width = 0.7, fill = "lightgrey") +
  geom_text(aes(label = n), vjust = -0.5) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            vjust = 1.5, color = "black", fontface = "bold") +
  
  scale_x_discrete(labels = c(
    "cds"        = "CDS",
    "utr5"       = "5' UTR",
    "utr3"       = "3' UTR",
    "introns"    = "Introns",
    "intergenic" = "Intergenic"
  )) +
  
  theme_minimal() +
  labs(title = "Genomic Distribution of Unique RNA Editing Sites",
       y = "Number of Sites",
       x = "") +
  
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black")
  )

print(p)


## Repetitive elements


tissue_names <- c("cortex", "hypothalamus", "hippocampus")
total_sites_count <- 15288  

rep_el_data_list <- lapply(tissue_names, function(tissue) {
  file <- paste0(tissue, '_repetitive_element_sites.bed')
  read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
})

rep_el_all_data <- do.call(rbind, rep_el_data_list)


rep_el_all_data$V5 <- trimws(rep_el_all_data$V5)
rep_el_unique_data <- rep_el_all_data[!duplicated(rep_el_all_data$V5), ]


rep_el_unique_data$Category <- ifelse(grepl("SINE", rep_el_unique_data$V12, ignore.case = TRUE), "SINE",
                               ifelse(grepl("LINE", rep_el_unique_data$V12, ignore.case = TRUE), "LINE", 
                                      "Other"))

rep_el_plot_data <- rep_el_unique_data %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / total_sites_count) * 100)


ggplot(rep_el_plot_data, aes(x = Category, y = Percentage)) + 
  geom_col(fill = "#56B4E9", width = 0.7) + 
  

  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            vjust = -0.5, 
            size = 5, 
            fontface = "bold") +
  
  scale_y_continuous(limits = c(0, 100)) + 
  labs(
    y = paste0('% of Total Unique Editing Sites (n=', total_sites_count,')'),
    x = NULL,
    title = 'Distribution of Editing Sites in Repetitive Elements'
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )


top_types_plot_data <- rep_el_unique_data %>%
  group_by(V11) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count)) %>%   # Sort highest to lowest
  slice_head(n = 6) %>%      # Take the top 7
  mutate(Percentage = (Count / total_sites_count) * 100)


ggplot(top_types_plot_data, aes(x = reorder(V11, -Percentage), y = Percentage)) + 
  geom_col(fill = "#E69F00", width = 0.7) +  
  
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            vjust = -0.5, 
            size = 5, 
            fontface = "bold") +
  
  scale_y_continuous(limits = c(0, max(top_types_plot_data$Percentage) + 5)) +
  labs(
    y = paste0('% of Total Unique Editing Sites (n=',total_sites_count, ')'),
    x = NULL,
    title = 'Top 7 Repetitive Elements in Editing Sites'
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1), 
    axis.title.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )


## editing rate distribution


all_rates <- as.vector(editing_matrix)

plot_editing_rate_distr <- tibble(rate = all_rates) %>%

  filter(!is.na(rate) & rate > 0) %>%
  mutate(
    bin = cut(rate, 
              breaks = seq(0, 1, by = 0.1), 
              include.lowest = FALSE, 
              right = TRUE,
              labels = c("(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]", 
                         "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]")
    )
  ) %>%
  dplyr::count(bin) %>%
  mutate(percentage = n / sum(n))


p <- ggplot(plot_editing_rate_distr, aes(x = bin, y = percentage)) +
  
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  
  # format y-axis as percentages
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  labs(
    title = "Distribution of RNA Editing Levels across editing events",
    x = "RNA editing levels",
    y = "Frequency"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p)


## Genes with highest number of editing sites, plot was not used for report

cds_sites <- read.table('genes_with_cds_editing_sites.tsv')
filtered_cds_sites <- cds_sites[grepl("synonymous|missense", cds_sites$V3), ]


top_15_list <- filtered_cds_sites %>%
  dplyr::count(V2) %>% 
  arrange(desc(n)) %>%
  slice_head(n = 15)

plot_data_cds_sites <- filtered_cds_sites %>%
  filter(V2 %in% top_15_list$V2) %>%
  dplyr::count(V2, V3) %>%
  mutate(V2= factor(V2, levels = rev(top_15_list$V2)))


ggplot(plot_data_cds_sites, aes(y = V2, x = n, fill = V3)) +
  geom_col() +
  labs(
    title = "Top 15 Genes by Editing Site Count",
    y = "Gene Name",
    x = "Number of Sites",
    fill = "Variant Type"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")


##bubble plot

##first extract the known Q607R editing site of GRIA2 to include it in the
#plot later

snps_across_tissues <- list(
  cortex = readRDS('all_SNPs_cortex_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hippocampus = readRDS('all_SNPs_hippocampus_rescue_with_add_filter_g_2_e_005_r_5.rds'),
  hypothalamus = readRDS('all_SNPs_hypothalamus_rescue_with_add_filter_g_2_e_005_r_5.rds')
)
  
long_snps <- snps_across_tissues %>%
  map(function(chrom_list) {
    
    chrom_list %>%
      map(function(rep_list) {
        
        bind_rows(rep_list, .id = "replicate_id")
        
      }) %>%
      
      list_rbind(names_to = "chromosome")
    
  }) %>%
  
  list_rbind(names_to = "tissue")

long_snps <- long_snps %>%
  mutate(unique_sample_id = paste(tissue, replicate_id, sep = "_"))


wide_snps <- long_snps %>%
  dplyr::select(unique_sample_id, site_id, editing_rate) %>%
  pivot_wider(
    names_from = site_id, 
    values_from = editing_rate,
    values_fill = 0 
  )

snp_matrix <- wide_snps %>% 
  column_to_rownames("unique_sample_id")


neuro_sites <- c('1:68530727-68530727:+', '1:68530740-68530740:+', '1:68564368-68564368:+',
                 '1:68564374-68564374:+', '1:68564381-68564381:+', '8:46206336-46206336:+',
                 '4:8670185-8670185:+','4:8670186-8670186:+', '4:8670331-8670331:+',
                 '5:65562460-65562460:-', '13:35884880-35884880:+', '13:35884888-35884888:+')


neuro_genes <- c('8:16915172(GRIA2)','8:46206336(GRIA2)', '5:65562460(KCNA1)',
                 '4:8670331(KCNQ3)', '4:8670186(KCNQ3)', '4:8670185(KCNQ3)',
                 '1:68564381(GRIK2)', '1:68564374(GRIK2)', '1:68564368(GRIK2)',
                 '1:68530740(GRIK2)', '1:68530727(GRIK2)', '13:35884888(CACNA1D)',
                 '13:35884880(CACNA1D)')

n_total <- nrow(editing_matrix)

# Create tissue vector
tissue_reps <- c(
  rep("Cortex", 48),
  rep("Hippocampus", 48),
  rep("Hypothalamus", 45)
)

neuro_edit_sub <- editing_matrix[, neuro_sites, drop = FALSE]


neuro_df_long <- as.data.frame(neuro_edit_sub) %>%
  mutate(
    sample_id = rownames(neuro_edit_sub),
    tissue = tissue_reps
  ) %>%
  pivot_longer(
    cols = all_of(neuro_sites),
    names_to = "site",
    values_to = "editing_rate"
  )


summary_df_neuro <- neuro_df_long %>%
  group_by(site, tissue) %>%
  summarise(
    n_samples = sum(editing_rate > 0, na.rm = TRUE),
    mean_editing_rate = mean(editing_rate[editing_rate > 0], na.rm = TRUE),
    .groups = "drop"
  )

gria2_site <- tibble(
  site = "8:16915172-16915172:+",
  tissue = c("Cortex", "Hippocampus", "Hypothalamus"),
  mean_editing_rate = c(1, 1, 1),
  n_samples = c(47, 48, 45)
)

summary_df_neuro <- bind_rows(summary_df_neuro, gria2_site)

summary_df_neuro$site <- factor(
  summary_df_neuro$site,
  levels = unique(summary_df_neuro$site)
)

df_key <- sub(":(\\d+)-.*", ":\\1", levels(summary_df_neuro$site))


label_key <- sub("\\(.*\\)", "", neuro_genes)


site_labels <- setNames(
  neuro_genes[match(df_key, label_key)],
  levels(summary_df_neuro$site)
)

ggplot(summary_df_neuro, aes(x = tissue, y = site)) +
  geom_point(aes(size = n_samples,
                 color = mean_editing_rate)) +
  scale_color_viridis_c(name = "Mean Editing Rate") +
  scale_size(range = c(1, 8),
             name = "Number of samples") +
  scale_y_discrete(labels = site_labels) +
  labs(y = "Site Coordinate") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),   # remove x-axis title
    axis.text.x = element_text(angle = 45, hjust = 1)
  )