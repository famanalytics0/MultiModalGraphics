#' MultiModalGraphics: An R Package for Graphical Integration of Multimodal Datasets
#'
#' The MultiModalGraphics package provides rich functionalities for combining different data types
#' into an integrative visual representation. This package facilitates the creation of intuitive visualizations
#' that merge multiple streams of data into a single coherent graphic.
#'
#' @name MultiModalGraphics-package
#' @docType package
#' @details
#' The MultiModalGraphics package offers the following features:
#' \itemize{
#'   \item InformativeHeatmap: encapsulating a CompexHeatmap object to provide additional functionality
#'   \item MultifeatureGrid: generating a heatmap grid, where each row presents a composite of information about a group of related features per row, rather than a single feature (unlike typically seen in conventional heatmaps).
#'   \item ClearScatterplot: generating customizable scatter plots that can encompass multiple datasets within a single visualization.
#'   \item generate_data: helper class to provide data to quickly test the 3 classes above.
#' }
#'
#' @section Getting Started:
#' To get started with the MultiModalGraphics package, you can install it from GitHub:
#' \preformatted{
#' # Install the development version from GitHub
#' # install.packages("devtools")
#' devtools::install_github("famanalytics0/MultiModalGraphics")
#' }
#'
#' @examples
#' # Example of creating an InformativeHeatmap
#' library(ComplexHeatmap)
#' library(MultiModalGraphics)
#' data <- matrix(rnorm(100), ncol = 10)
#' heatmap <- InformativeHeatmap(data, pch_val = 20, unit_val = 2,
#'                                 significant_color = "red",
#'                                 trending_color = "blue",
#'                                 significant_pvalue = 0.05,
#'                                 trending_pvalue = 0.1)
#' draw(getHeatmapObject(heatmap))
#'
#' # Example of creating a ClearScatterplot
#' plotdata <- get_clear_scatterplot_df()
#' scatterplotObject <- ClearScatterplot(
#'   data = plotdata,
#'   logFoldChange = "log2fc",
#'   timePointColumn = "timePoint",
#'   timePointLevels = c("T10R1", "T5R1")
#' )
#' scattered_plot <- createPlot(
#'   scatterplotObject,
#'   color1 = "cornflowerblue",
#'   color2 = "grey",
#'   color3 = "indianred",
#'   highLog2fc = 0.585,
#'   lowLog2fc = -0.585,
#'   negLog10pValue = 1.301,
#'   expressionDirection = "regulation",
#'   negativeLogPValue = "negLog10p",
#'   timeVariable = "reg_time_org",
#'   xAxis = "organ",
#'   yAxis = "timePoint"
#' )
#' print(scattered_plot)
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ggplot2 ggplot geom_point geom_text geom_jitter facet_grid theme scale_fill_hue scale_color_manual labs element_text element_rect
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr %>% group_by summarize
#' @import methods
NULL
