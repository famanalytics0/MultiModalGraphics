# Defines the S4 class 'ClearScatterplot' to represent scattered plot with specific metadata and plot configurations
#' ClearScatterplot: An S4 class for scatterplot visualization
#'
#' This class represents a scatterplot with specific metadata and plot configurations.
#'
#' @slot metadata_filepath Character vector: Path to metadata file.
#' @slot metadata Data.frame: Metadata for the scatterplot.
#' @slot plot ANY: The ggplot2 plot object (initially NULL).
#' @name ClearScatterplot
#' @docType class
#' @importFrom methods setClass
#' @export
setClass("ClearScatterplot",
         representation(
           metadata_filepath = "character",
           metadata = "data.frame",
           plot = "ANY"
         )
)

#' @title Initialize a ClearScatterplot Object
#' @description Constructor method for the `ClearScatterplot` class.
#'   Initializes a `ClearScatterplot` object with specified metadata and plot configurations.
#' @param .Object The `ClearScatterplot` object to be initialized.
#' @param metadata_filepath Character string specifying the path to the metadata file.
#'   The metadata should contain information necessary for creating scatter plots.
#' @param timePointLevels Character vector specifying the levels for time points.
#'   Defaults to c("TP-1", "TP-2", "TP-3", "TP-4").
#' @param pValueColumn Character string specifying the name of the column in the metadata
#'   that contains p-values. Default is "p".
#' @param qValueColumn Character string specifying the name of the column in the metadata
#'   that contains q-values. Default is "q".
#' @param expressionColumnName Character string specifying the name of the column in the metadata
#'   that contains expression values. Default is "log2fc".
#' @param highLog2fc Numeric value specifying the threshold for high log2 fold change.
#'   Default is 0.585.
#' @param lowLog2fc Numeric value specifying the threshold for low log2 fold change.
#'   Default is -0.585.
#' @param negLog10pValue Numeric value specifying the threshold for negative log10 p-values.
#'   Default is 1.301.
#' @return An initialized `ClearScatterplot` object with metadata loaded and ready for plot creation.
#' @export
#' @examples
#' \dontrun{
#'   metadata_path <- system.file("extdata", "example_metadata.csv", package = "YourPackage")
#'   csPlot <- new("ClearScatterplot", metadata_filepath = metadata_path)
#' }
setMethod("initialize",
          signature(.Object = "ClearScatterplot"),
          function(.Object, metadata_filepath, timePointLevels = c("TP-1", "TP-2", "TP-3", "TP-4"),
                   pValueColumn = "p", qValueColumn = "q", expressionColumnName = "log2fc",
                   highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301) {
            # Initialize a 'ClearScatterplot' object with the given metadata and default or specified parameters
            .Object@metadata_filepath <- metadata_filepath # Path to metadata file
            .Object@metadata <- process_metadata(metadata_filepath, pValueColumn = pValueColumn,
                                                 qValueColumn = qValueColumn, expressionColumnName = expressionColumnName,
                                                 highLog2fc = highLog2fc, lowLog2fc = lowLog2fc, negLog10pValue =  negLog10pValue)
            .Object@plot <- NULL # Placeholder for the plot, to be created later
            return(.Object) # Return the initialized object
          }
)

# Define a generic method 'createPlot' to create plots for objects
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' @title Create and Assign a Plot to a ClearScatterplot Object
#' @description Generates a scatter plot based on the stored metadata and aesthetic parameters provided. This plot is then assigned to the `plot` slot of the `ClearScatterplot` object.
#' @param object A `ClearScatterplot` object for which the plot will be created.
#' @param color1 Character string specifying the color for high significance points. Default is "cornflowerblue".
#' @param color2 Character string specifying the color for points that do not meet the high or low significance thresholds. Default is "grey".
#' @param color3 Character string specifying the color for low significance points. Default is "indianred".
#' @param highLog2fc Numeric threshold for considering an expression value significantly high. Default is 0.585.
#' @param lowLog2fc Numeric threshold for considering an expression value significantly low. Default is -0.585.
#' @param expressionDirection Character string indicating the direction of regulation to be highlighted. Possible values are "up", "down", or "regulation" for both. Default is "regulation".
#' @param timeVariable Character string specifying the column in the metadata that represents the time variable. This parameter can be used to facet the plot by time points. Default is "reg_time_org".
#' @return The `ClearScatterplot` object with its `plot` slot updated to include the newly created plot.
#' @export
#' @examples
#' \dontrun{
#'   metadata_path <- system.file("extdata", "example_metadata.csv", package = "YourPackage")
#'   csPlot <- new("ClearScatterplot", metadata_filepath = metadata_path)
#'   csPlot <- createPlot(csPlot)
#'   # Now, csPlot contains a ggplot object in its plot slot that can be displayed using show(csPlot)
#' }
setMethod("createPlot",
          signature(object = "ClearScatterplot"),
          function(object, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                   highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                   timeVariable="reg_time_org") {
            # Create a plot based on the metadata and specified aesthetic parameters
            if (is.null(object@plot)) {
              object@plot <- create_plot(object@metadata, color1 = color1, color2 = color2, color3=color3,
                                         highLog2fc = highLog2fc, lowLog2fc = lowLog2fc,
                                         expressionDirection = expressionDirection,
                                         timeVariable = timeVariable)
            }
            return(object) # Return the object with its plot
          }
)

# Function to create the actual plot based on processed metadata and aesthetic choices
create_plot <- function(metadata, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                        highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                        timeVariable="reg_time_org") {
  xAxis = "organ"
  yAxis = "timePoint"
  colorFlag = "color_flag"
  facetFormula = paste0(xAxis, "~", yAxis)

  gp_obj <- ggplot(data = metadata, ggplot2::aes(
    x = log2fc, y = neglog10p, color = as.factor(metadata[[colorFlag]])
  )) +
    ggplot2::geom_point(alpha = 0.5, size = 1.75) +
    ggplot2::theme(legend.position = "none") + ggplot2::geom_jitter() +
    ggplot2::scale_color_manual(values = c(color1, color2, color3)) +
    ggplot2::labs(
      x = expression("log2 (fold change)"), y = expression("-log10 (p-value)")
    ) + ggplot2::facet_grid(facetFormula, space = "free") + ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12)
    )

  gp_obj2 = gp_obj + ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),axis.title.x = ggplot2::element_text(size = 12,
                                                                                      face = "bold"),title = ggplot2::element_text(size = 12, face = "bold"),
    strip.text.x = ggplot2::element_text(size = 12, face = "bold",colour = "black", angle = 0)
  )
  gp_obj3 = gp_obj2 + ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),axis.title.y = ggplot2::element_text(size = 12,
                                                                                      face = "bold"),title = ggplot2::element_text(size = 12, face = "bold"),
    strip.text.y = ggplot2::element_text( size = 12, face = "bold", colour = "black", angle = 90)
  )
  p <- gp_obj3 + ggplot2::theme(legend.position = "right")

  p1 <-
    p + ggplot2::theme(strip.background = ggplot2::element_rect(
      fill = "white", color = color2, linewidth = 1
    ))

  p2 <- p1 + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::scale_color_manual(name = "log2 (fold change)", values = c(color1, color2, color3),
                       labels = c(paste("<", lowLog2fc), paste(lowLog2fc,"to", highLog2fc), paste(">", highLog2fc))
    ) + ggplot2::scale_fill_hue()

  metadata2 <- metadata[metadata[[colorFlag]] == 1, ]
  metadata3 <- metadata[metadata[[colorFlag]] == -1, ]
  metadata4 <- rbind(metadata2, metadata3)

  metadata5 <- metadata4[metadata4[[expressionDirection]] == "up", ]
  metadata6 <- metadata4[metadata4[[expressionDirection]] == "down", ]

  # meta.cor1 <- plyr::ddply(.data = metadata5,
  #                          .(organ, timePoint),
  #                          .fun = function(x) {
  #                            dplyr::summarize(x, n1 = paste(length(x[[timeVariable]])))
  #                          }
  # )

  # meta.cor1 <- dplyr::summarize(
  #   dplyr::group_by(metadata5, metadata5[[xAxis]], metadata5[[yAxis]]),
  #   n1 = paste(length(metadata[[timeVariable]])),
  #   .groups = 'drop'
  # )

  meta.cor1 <- metadata5 %>%
    dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) %>%
    dplyr::summarize(n1 = paste(length(.data[[timeVariable]])), .groups = 'drop')

  p3 <- p2 + ggplot2::geom_text(data = meta.cor1, ggplot2::aes(x = 4, y = 9, label = n1),
                                colour = color3, inherit.aes = FALSE, parse = TRUE
  )

  # meta.cor2 <- plyr::ddply(.data = metadata6,
  #                          .(organ, timePoint),
  #                          .fun = function(x) {
  #                            dplyr::summarize(x, n2 = paste(length(x[[timeVariable]])))
  #                          }
  # )

  meta.cor2 <- metadata6 %>%
    dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) %>%
    dplyr::summarize(n2 = paste(length(.data[[timeVariable]])), .groups = 'drop')


  # meta.cor2 <- dplyr::summarize(
  #   dplyr::group_by(metadata6, metadata6[[xAxis]], metadata6[[yAxis]]),
  #   n2 = paste(length(metadata[[timeVariable]])),
  #   .groups = 'drop'
  # )

  p4 <- p3 + ggplot2::geom_text(data = meta.cor2, ggplot2::aes(x = -4, y = 9, label = n2),
                                colour = color1, inherit.aes = FALSE, parse = FALSE
  )

  # Return the created plot
  return(p4 + ggplot2::theme(legend.position = "bottom"))
}


#' @title Display the Scatterplot Stored in a ClearScatterplot Object
#' @description This method invokes the plotting functionality to generate and display the scatterplot based on the metadata and configuration stored within the `ClearScatterplot` object. If the plot has not been created yet, it first calls the `createPlot` method to generate the plot and then displays it.
#' @param object A `ClearScatterplot` object containing the necessary metadata and configuration for generating the scatterplot.
#' @return Invisible `ClearScatterplot` object. The primary purpose of this method is the side effect of displaying the plot, not returning a value. However, the object itself is returned invisibly for potential further manipulation or inspection.
#' @export
#' @examples
#' \dontrun{
#'   metadata_path <- system.file("extdata", "example_metadata.csv", package = "YourPackage")
#'   csPlot <- new("ClearScatterplot", metadata_filepath = metadata_path)
#'   # Assuming createPlot method has been properly defined and exported
#'   show(csPlot)
#'   # This will display the scatterplot based on the metadata
#' }
setMethod("show",
          signature(object = "ClearScatterplot"),
          function(object) {
            object <- createPlot(object)
            print(object@plot)
          }
)

# Function to process metadata for 'ClearScatterplot' objects
process_metadata <- function(filepath, timePointLevels, pValueColumn = "p", qValueColumn = "q",
                             expressionColumnName = "log2fc", highLog2fc = 0.585, lowLog2fc = -0.585,
                             negLog10pValue =  1.301) {
  metadata <- read.csv(filepath)
  metadata$neglog10p = -log10(metadata[[pValueColumn]])
  #metadata$neglog10q = -log10(metadata[[qValueColumn]])

  #this for p < 0.05 and FC > 1.5 (linear fold change)
  metadata$color_flag <-
    ifelse(
      metadata$log2fc > highLog2fc &
        metadata$neglog10p > negLog10pValue, 1, ifelse(metadata$log2fc < lowLog2fc &
                                                         metadata$neglog10p > negLog10pValue,-1, 0)
    )

  metadata$timePoint <- factor(
    metadata$timePoint, levels = timePointLevels
  )
  return(metadata)
}
