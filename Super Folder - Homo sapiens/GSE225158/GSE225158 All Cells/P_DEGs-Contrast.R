library(ggplot2)
library(dplyr)
library(stringr)

# List all DESeq2 result files
contrast_files <- list.files(
  path = '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/',
  pattern = '^deseq2_results_.*\\.csv$',
  full.names = TRUE
)

# Initialize list to store DEG counts
deg_counts_list <- list()

# Loop through each file, count DEGs, and store results
for (file in contrast_files) {
  res <- read.csv(file)
  # Define DEGs: padj < 0.05 and abs(log2FoldChange) > 1
  degs <- res %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  # Extract contrast name from filename
  contrast_name <- str_replace(basename(file), "^deseq2_results_|\\.csv$", "")
  deg_counts_list[[contrast_name]] <- nrow(degs)
}

# Create summary data frame
deg_counts_df <- data.frame(
  contrast = names(deg_counts_list),
  n_DEGs = as.integer(unlist(deg_counts_list))
)

# Barplot
ggplot(deg_counts_df, aes(x = contrast, y = n_DEGs, fill = contrast)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of DEGs per Contrast", x = "Contrast", y = "Number of DEGs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Save the plot as a PNG file
ggsave("DEGs_per_contrast.png", width = 8, height = 6, dpi = 300)