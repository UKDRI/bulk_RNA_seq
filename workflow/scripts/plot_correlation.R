library(ggplot2)
library(reshape2)

# Load FPKM data
fpkm <- read.csv(snakemake@input[[1]], row.names = 1)

# Compute correlation matrix
cor_matrix <- cor(fpkm, method = "pearson")
melted_cor <- melt(cor_matrix)

# Plot heatmap
p <- ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  labs(title = "Pearson Correlation Between Samples", x = "", y = "", fill = "R^2") +
  theme_minimal()

# Save plot
ggsave(snakemake@output[[1]], plot = p, width = 8, height = 6)
