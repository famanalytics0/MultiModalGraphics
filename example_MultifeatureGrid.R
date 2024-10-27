library(MultiModalGraphics)
# Example usage:
filepath <- "data_MultifeatureGrid/twoD_graphics.csv"
data <- read.csv(filepath)

# Create a MultifeatureGrid object
mg <- MultifeatureGrid(data,breaks = seq(-3, 3, 1))

plot_heatmap(mg,  independantVariable = "timePoint")


# Plot the heatmap

heatmap_obj <- MultifeatureGrid(data, title = "Trauma Exposure and Post-Trauma Days",
                                x_label = "Brain Region",
                                y_label = "Neuronal Signaling, Synaptic Plasticity and Neurogenesis",
                                logpval_label = "-log10(p-value)",
                                zscore_label= "Activation z-score",
                                numitems_label ="Number of Genes",
                                color_palette = "PuOr", breaks = seq(-3, 3, 1))

# Plot the heatmap
plot_heatmap(heatmap_obj, pValueColumn = "p", lowColor = "yellow", highColor = "red",
             borderColor="grey60", columnForNumber = "number_of_genes",  independantVariable = "timePoint")
