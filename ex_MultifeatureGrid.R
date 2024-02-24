# Example usage:
# Create a MultifeatureGrid object

# setwd("/cloud/project/rstudio")
# Load or source the MultifeatureGrid.R file
data <- data.frame(tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
                   signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
                   Activation_z_score = runif(8, -2, 2),
                   p = runif(8, 0, 0.05),
                   number_of_genes = sample(1:100, 8))
data
mg <- MultifeatureGrid(data)
plot_heatmap(mg,  independantVariable = "tissue")

# Plot the heatmap

