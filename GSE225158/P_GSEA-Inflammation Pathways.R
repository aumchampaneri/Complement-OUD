# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

# This dot plot visualizes GSEA pathway enrichment results:
#
# Y-axis: Pathway names (each row is a pathway).
# X-axis: Normalized Enrichment Score (NES). Higher absolute NES indicates stronger enrichment (positive = upregulated, negative = downregulated).
# Dot color: Reflects NES value (blue = negative, red = positive, white = neutral).
# Dot size: Represents statistical significance (-log10(padj)); larger dots are more significant.
# Interpretation:
# Look for large, colored dots far from zero on the x-axis—these are the most strongly and significantly enriched pathways in your data. Pathways with large red dots (high positive NES) are upregulated; large blue dots (high negative NES) are downregulated.


# Read and preprocess the data
df <- read.csv('/GSE225158/GSEA outputs/gsea_inflammatory_complement_pathways.csv',
               header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("Name", "Description", "Size", "ES", "NES", "pval", "padj", "qval", "Rank", "Tags", "Genes")
df$Name <- str_replace_all(df$Name, '"', '')
df$NES <- as.numeric(df$NES)
df$padj <- as.numeric(df$padj)

# Order by NES for better visualization
df <- df %>% arrange(NES)

# Create the dot plot
p <- ggplot(df, aes(x = NES, y = reorder(Name, NES), color = NES, size = -log10(padj))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(padj)",
    title = "GSEA Pathway Enrichment"
  ) +
  theme(axis.text.y = element_text(size = 8))

# Save the plot with increased width and left margin
ggsave(
  '/Users/aumchampaneri/PycharmProjects/Complement-OUD/gsea_inflammatory_complement_pathways_dotplot.png',
  p + theme(plot.margin = margin(5, 5, 5, 20)), # increase left margin
  width = 20, height = 8 # increase width
)