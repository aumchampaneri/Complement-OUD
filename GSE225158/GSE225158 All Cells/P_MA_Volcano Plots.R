# R

library(ggplot2)

contrast_names <- c(
  "F_OUD_vs_F_None",
  "M_OUD_vs_M_None",
  "F_OUD_vs_M_OUD",
  "F_None_vs_M_None"
)

for (contrast in contrast_names) {
  res_path <- paste0('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/deseq2_results_', contrast, '.csv')
  res <- read.csv(res_path)

  # ... inside your loop, before the MA plot:
res_ma <- res[res$baseMean > 0, ]
png(paste0('MA_plot_', contrast, '.png'), width=800, height=600)
with(res_ma, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log='x',
     main=paste('MA Plot:', contrast), xlab='Mean Expression', ylab='log2 Fold Change'))
abline(h=0, col='red')
dev.off()

  # MA Plot
  png(paste0('MA_plot_', contrast, '.png'), width=800, height=600)
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log='x',
       main=paste('MA Plot:', contrast), xlab='Mean Expression', ylab='log2 Fold Change'))
  abline(h=0, col='red')
  dev.off()

  # Volcano Plot
  res$padj[is.na(res$padj)] <- 1
  res$significant <- res$padj < 0.05
  png(paste0('Volcano_plot_', contrast, '.png'), width=800, height=600)
  p <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c('grey','red')) +
    theme_minimal() +
    labs(title=paste('Volcano Plot:', contrast), x='log2 Fold Change', y='-log10 Adjusted p-value')
  print(p)  # <-- This line is required
  dev.off()
}