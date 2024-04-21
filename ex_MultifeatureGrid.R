library(MultiModalGraphics)
# Example usage:
# Create a MultifeatureGrid object
data <- get_multifeature_grid_df()
mg <- MultifeatureGrid(data)
plot_heatmap(mg,  independantVariable = "tissue")

# Plot the heatmap

