# MultiModalGraphics Classes

# 1. ClearScatterplot

## Introduction

This  example demonstrates how to use the `ClearScatterplot` class to create 
enhanced scatterplots. The `ClearScatterplot` class provides functionalities 
for creating scatterplots with enhanced visualization features like coloring 
based on significance levels and plotting against different variables.

## Example: Creating a Scatterplot

First, load the necessary data and create an instance of `ClearScatterplot`.

```{r load-data}
plotdata <- get_clear_scatterplot_df()

scatterplotObject <- ClearScatterplot(
  data = plotdata,
  logFoldChange = "log2fc",
  timePointColumn = "timePoint",
  timePointLevels = c("T10R1", "T5R1")
)
```

Next, create the scatterplot using the `createPlot` method, specifying various 
parameters such as color and thresholds.

```{r create-plot}
scattered_plot <-
  createPlot(
    scatterplotObject,
    color1 = "cornflowerblue",
    color2 = "grey",
    color3 = "indianred",
    highLog2fc = 0.585,
    lowLog2fc = -0.585,
    negLog10pValue = 1.301,
    expressionDirection = "regulation",
    negativeLogPValue = "negLog10p",
    timeVariable = "reg_time_org",
    xAxis = "organ",
    yAxis = "timePoint"
  )
```

Finally, display the plot. This will call the 'show' method to render the plot.

```{r display-plot}
scattered_plot  # Display the plot
```

## Session info
```{r, echo=FALSE}
sessionInfo()
```

# 2. InformativeHeatmap

This section demonstrates how to create and visualize an informative heatmap 
using the `MultiModalGraphics` package which utilizes the `ComplexHeatmap` 
package for enhanced visualizations.

## Basic Usage

Below is an example of how to create an informative heatmap with data 
representing genes, their value groups, and significance levels.

### Data Preparation

First, we load the necessary data. 

```{r load-data}
informative_heatmap <- get_informative_heatmap_df()
informative_heatmap_matrix <- as.matrix(informative_heatmap)
group_val <- informative_heatmap_matrix[, 1:3]
p_val <- informative_heatmap_matrix[, 4:6]
```

### Creating Heatmap

Create an `InformativeHeatmap` object with custom settings for visual 
representation.

```{r create-heatmap}
htmp_plus <- InformativeHeatmap(group_val,
                                unit_val = 7,
                                pch_val = 16,
                                significant_color = "black",
                                trending_color = "turquoise",
                                significant_pvalue = 0.05,
                                trending_pvalue = 0.1,
                                significance_level = p_val,
                                row_title = "Genes",
                                column_title = "Value and Significance",
                                cluster_rows = TRUE,
                                show_row_names = TRUE,
                                row_names_side = "left",
                                column_names_rot = 45,
                                row_dend_reorder = TRUE,
                                rect_gp = gpar(col = "white", lwd = 2))
```

### Drawing Heatmap

Finally, draw the heatmap using the ComplexHeatmap function to visualize the 
data.

```{r draw-heatmap}
draw(getHeatmapObject(htmp_plus))
```

## Conclusion

This vignette provided a simple example of how to use the `MultiModalGraphics` 
package to create a detailed heatmap for data analysis. For more advanced 
features and customization, refer to the `ComplexHeatmap` documentation.

## Session info
```{r, echo=FALSE}
sessionInfo()
```


# 3. MultifeatureGrid

This example demonstrates how to use the `MultifeatureGrid` class from the 
`MultiModalGraphics` package to create comprehensive heatmaps integrating 
various data features such as z-scores, p-values, and counts.

## Example Data Preparation

First, we create a sample dataset representing biological data across different 
tissues and signaling pathways, with associated p-values and activation 
z-scores.

```{r example-data}
data <- get_multifeature_grid_df()

```

## Creating a MultifeatureGrid Object

We initialize the `MultifeatureGrid` object with the data prepared above.

```{r create-grid}
mg <- MultifeatureGrid(data)
```

## Plotting the Heatmap

We then plot the heatmap, specifying 'tissue' as the independent variable for 
faceting.

```{r plot-heatmap}
plot_heatmap(mg, independantVariable = "tissue")
```

This plot provides a visual summary of the signaling activity and the 
statistical significance across different tissues, utilizing a color gradient 
to represent activation z-scores and the size of points to indicate the number 
of genes involved.

## Session info
```{r, echo=FALSE}
sessionInfo()
