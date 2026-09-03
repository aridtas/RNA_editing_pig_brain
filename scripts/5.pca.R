library(tidyverse)
library(ggplot2)


constant_columns <- apply(editing_matrix, 2, var, na.rm = TRUE) == 0
pca_matrix_clean <- editing_matrix[ , !constant_columns]



pca_res <- prcomp(editing_matrix, center = TRUE, scale. = TRUE)



pca_plot_data <- as.data.frame(pca_res$x) %>%
  rownames_to_column("unique_sample_id") %>%
  # Extract back the Tissue info from the ID string for coloring
  separate(unique_sample_id, into = c("Tissue", "Replicate"), sep = "_", 
                                      remove = FALSE)

# Proportion of variance explained
var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of RNA Editing Profiles",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "right")




# ggsave("PCA_editing_sites.png", plot = p, width = 8, height = 6)

