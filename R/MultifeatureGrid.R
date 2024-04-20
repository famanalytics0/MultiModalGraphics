#' MultifeatureGrid Class for 2D Heatmap Visualization
#'
#' This class represents a 2D heatmap with customizable configurations for plotting
#' biological data, integrating features like significance levels and z-scores.
#'
#' @slot data Data frame containing heatmap data.
#' @slot title Character vector specifying the heatmap title.
#' @name MultifeatureGrid
#' @docType class
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradientn geom_point scale_color_gradient scale_size labs facet_grid theme_bw theme element_text element_rect
#' @importFrom RColorBrewer brewer.pal
#' @export MultifeatureGrid
#' @export plot_heatmap
#' @examples
#' data <- data.frame(tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
#'                    signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
#'                    Activation_z_score = runif(8, -2, 2),
#'                    p = runif(8, 0, 0.05),
#'                    number_of_genes = sample(1:100, 8))
#' mg <- MultifeatureGrid(data)
#' plot_heatmap(mg)
setClass(
  "MultifeatureGrid",
  slots = list(
    data = "data.frame",          # Data frame containing the heatmap data
    title = "character",          # Title of the heatmap
    x_label = "character",        # Label for the X-axis
    y_label = "character",        # Label for the Y-axis
    logpval_label = "character",  # Label for the Log(p-value) legend
    zscore_label = "character",   # Label for the z-score legend
    numitems_label = "character", # Label for the number of items legend
    color_palette = "character",  # Name of the color palette to use
    breaks = "numeric"            # Numeric vector specifying the breakpoints for color mapping
  )
)

#' Constructor for MultifeatureGrid Class
#'
#' This function initializes a `MultifeatureGrid` object with data and configuration
#' for creating a multi-feature grid heatmap. It allows customization of the heatmap's
#' title, axis labels, color palette, and more.
#'
#' @param data A `data.frame` containing the data to be visualized in the heatmap.
#'   The data frame should have columns that match the expected metadata fields
#'   (e.g., significance values, z-scores).
#' @param title The title of the heatmap. Default is "Heatmap".
#' @param x_label The label for the x-axis. Default is "X Label".
#' @param y_label The label for the y-axis. Default is "Y Label".
#' @param logpval_label The label for the legend describing the -log10(p-value).
#'   Default is "-log10(p-value)".
#' @param zscore_label The label for the legend describing the activation z-score.
#'   Default is "Activation z-score".
#' @param numitems_label The label for the legend describing the number of genes
#'   or items. Default is "Number of Genes".
#' @param color_palette The name of the color palette to use for the heatmap.
#'   Default is "RdYlBu". This should be a valid palette name recognized by RColorBrewer.
#' @param breaks A numeric vector specifying the breakpoints for color mapping
#'   across the heatmap's values. Default is `seq(-1, 1, 0.5)`.
#' @return An object of class `MultifeatureGrid` initialized with the specified
#'   parameters.
#' @export
#' @examples
#' \dontrun{
#'   data <- data.frame(
#'     significance = runif(100),
#'     z_score = rnorm(100),
#'     num_genes = sample(1:100, 100, replace = TRUE)
#'   )
#'   mg <- MultifeatureGrid(data)
#'   # Further usage of mg...
#' }
MultifeatureGrid <- function(data, title = "Heatmap", x_label = "X Label",
                             y_label = "Y Label", logpval_label = "-log10(p-value)",
                             zscore_label= "Activation z-score",
                             numitems_label ="Number of Genes",
                             color_palette = "RdYlBu", breaks = seq(-1, 1, 0.5)) {
  # Initializes a MultifeatureGrid object with provided or default parameters
  new("MultifeatureGrid", data = data, title = title, x_label = x_label, y_label = y_label,
      logpval_label = logpval_label, zscore_label = zscore_label,
      numitems_label = numitems_label, color_palette = color_palette, breaks = breaks)
}

# Method for plotting the heatmap for MultifeatureGrid objects
setGeneric("plot_heatmap", function(object, ...) standardGeneric("plot_heatmap"))

#' Plot Heatmap for MultifeatureGrid Objects
#'
#' This method generates and plots a heatmap for `MultifeatureGrid` objects. It allows
#' customization of the heatmap's appearance based on several parameters.
#'
#' @param object An object of class `MultifeatureGrid`. This object contains
#'   the data and configuration for the heatmap to be plotted.
#' @param pValueColumn A character string specifying the name of the column in the
#'   data frame that contains p-values. These values are used to calculate
#'   and plot the negative log10 of the p-values. Default is "p".
#' @param lowColor A character string specifying the color to use for the low
#'   end of the color gradient for p-values. Default is "yellow".
#' @param highColor A character string specifying the color to use for the high
#'   end of the color gradient for p-values. Default is "red".
#' @param borderColor A character string specifying the color to use for the border
#'   of each tile in the heatmap. Default is "grey60".
#' @param columnForNumber A character string specifying the name of the column in the
#'   data frame that contains the number of items (e.g., genes) for each feature.
#'   This value is used to determine the size of the points plotted on the heatmap.
#'   Default is "number_of_genes".
#' @param independantVariable A character string specifying the name of the column used
#'   as an independent variable for faceting the plot into multiple panels, one for each
#'   level of the independent variable. Default is "timePoint".
#'
#' @return A ggplot object representing the heatmap, which is also printed to the current
#'   plotting device.
#' @export
#' @examples
#' data <- data.frame(tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
#'                    signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
#'                    Activation_z_score = runif(8, -2, 2),
#'                    p = runif(8, 0, 0.05),
#'                    number_of_genes = sample(1:100, 8),
#'                    timePoint = rep(c("Time1", "Time2"), 4))
#' mg <- MultifeatureGrid(data)
#' plot_heatmap(mg)
setMethod("plot_heatmap", signature(object = "MultifeatureGrid"),
          function(object, pValueColumn = "p", lowColor = "yellow", highColor = "red",
                   borderColor="grey60", columnForNumber = "number_of_genes",
                   independantVariable = "timePoint") {
            # Extracts parameters from the object
            data <- object@data
            title <- object@title
            x_label <- object@x_label
            y_label <- object@y_label
            logpval_label <- object@logpval_label
            zscore_label <- object@zscore_label
            numitems_label <- object@numitems_label
            color_palette <- object@color_palette
            breaks <- object@breaks

            # Calculate negative log10 of p-values in the dataset
            data$neglog10p <- -log10(data[[pValueColumn]])

            # Define color palette
            color <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n <- 7, name <- color_palette)))(100)

            # Create the base ggplot object
            base_plot <- ggplot2::ggplot(data, ggplot2::aes(y <- signaling, x = tissue)) +
              ggplot2::geom_tile(ggplot2::aes(fill = Activation_z_score), colour = borderColor) +
              ggplot2::scale_fill_gradientn(colours = color, breaks = breaks, labels = scales::comma_format()) +
              ggplot2::geom_point(ggplot2::aes(colour = neglog10p, size = .data[[columnForNumber]])) +
              ggplot2::scale_color_gradient(low = lowColor, high = highColor) +
              ggplot2::scale_size(range = c(1, 10)) +
              ggplot2::labs(
                x <- x_label,
                y <- y_label,
                title <- title,
                fill <- zscore_label,
                colour <- logpval_label,
                size <- numitems_label
              ) +
              ggplot2::facet_grid(rows <- . ~ data[[independantVariable]], scales = "free_x", space = "free") +
              ggplot2::theme_bw() +
              ggplot2::theme(
                axis.text.x <- ggplot2::element_text(angle <- 90, hjust <- 1),
                panel.border <- ggplot2::element_rect(fill <- NA, colour <- "grey80", linewidth = 0.6),
                axis.text <- ggplot2::element_text(size <- 14, face <- "bold"),
                axis.title <- ggplot2::element_text(size <- 18, face <- "bold"),
                title <- ggplot2::element_text(size <- 18),
                strip.text.x <- ggplot2::element_text(size <- 14, face <- "bold", colour <- "black", angle <- 0)
              )

            # Display the plot
            print(base_plot)
          })
