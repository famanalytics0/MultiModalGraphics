#' ClearScatterplot Class
#'
#' A class for creating and managing scatterplot objects with enhanced visualization features.
#' @slot data A data frame containing the data to be plotted.
#' @slot plot An object storing the ggplot representation of the data.
#' @export
setClass("ClearScatterplot",
         slots = c(
           data = "data.frame",  # Slot for storing the data frame
           plot = "ANY"  # Slot for storing the ggplot object
         ))

setGeneric("ClearScatterplot", function(data, ...) {
  standardGeneric("ClearScatterplot")
})
#' ClearScatterplot Constructor
#'
#' Creates an instance of the ClearScatterplot class.
#' @param data A data frame containing the plot data.
#' @param logFoldChange The name of the column containing expression values.
#' @param negativeLogPValue The name of the column containing the negative log pValues
#' @param highLog2fc Threshold for high log2 fold change values.
#' @param lowLog2fc Threshold for low log2 fold change values.
#' @param negLog10pValue Threshold for -log10 p-value.
#' @param timePointColumn The name of the column containing time point information.
#' @param timePointLevels The levels for the time point column, if any.
#' @return An object of class ClearScatterplot.
#' @export
#' @examples
#' data <- data.frame(timePoint = c("T0", "T1"), p = runif(10), log2fc = runif(10, -2, 2))
#' scatterplot <- ClearScatterplot(data = data, pValueColumn = "p", expressionColumnName = "log2fc")
ClearScatterplot <- function(data, logFoldChange = "log2fc",
                             negativeLogPValue = "negLog10p",
                             highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue = 1.301,
                             timePointColumn = "timePoint", timePointLevels = NULL,
                             ...) {
  if (!is.data.frame(data)) {
    stop("Data must be a data frame.")
  }

  # Data preprocessing
  data$color_flag <- ifelse(
    data[[logFoldChange]] > highLog2fc & data[[negativeLogPValue]] > negLog10pValue, 1,
    ifelse(data[[logFoldChange]] < lowLog2fc & data[[negativeLogPValue]] > negLog10pValue, -1, 0)
  )

  # Factorize the timePoint column based on specified levels if provided
  if (!is.null(timePointLevels) && !is.null(data[[timePointColumn]])) {
    data[[timePointColumn]] <- factor(data[[timePointColumn]], levels = timePointLevels)
  }

  # Create the object
  obj <- new("ClearScatterplot", data = data)

  return(obj)
}


# Define a generic method 'createPlot' to create plots for objects
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' Create Plot Method for ClearScatterplot
#'
#' Creates a plot based on the ClearScatterplot object data.
#' @param object A ClearScatterplot object.
#' @param color1 Color for one category of data points.
#' @param color2 Color for another category of data points.
#' @param color3 Color for a third category of data points.
#' @param highLog2fc Threshold for high log2 fold change values.
#' @param lowLog2fc Threshold for low log2 fold change values.
#' @param expressionDirection Direction of gene expression.
#' @param negativeLogPValue The name of the column containing the negative log pValues
#' @param timeVariable The variable representing time.
#' @return The ClearScatterplot object with the plot updated.
#' @export
setMethod("createPlot",
          signature(object = "ClearScatterplot"),
          function(object, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                   highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue = 1.301,
                   expressionDirection = "regulation", negativeLogPValue="negLog10p",
                   timeVariable="reg_time_org",
                   xAxis = xAxis, yAxis = yAxis) {
            # Create a plot based on the data and specified aesthetic parameters
            if (is.null(object@plot)) {
              object@plot <- create_plot(object@data, color1 = color1, color2 = color2, color3=color3,
                                         highLog2fc = highLog2fc, lowLog2fc = lowLog2fc,
                                         expressionDirection = expressionDirection,
                                         negativeLogPValue = negativeLogPValue,
                                         timeVariable = timeVariable, xAxis = xAxis,
                                         yAxis = yAxis)
            }
            return(object) # Return the object with its plot
          }
)

# Function to create the actual plot based on processed data and aesthetic choices
create_plot <- function(data, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                        highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                        negativeLogPValue = "negLog10p", logFoldChange = "log2fc",
                        timeVariable="reg_time_org", xAxis = "organ", yAxis = "timePoint") {

  colorFlag = "color_flag"
  facetFormula = paste0(xAxis, "~", yAxis)

  gp_obj <- ggplot(data = data, ggplot2::aes(
    x = .data[[logFoldChange]], y = .data[[negativeLogPValue]], color = as.factor(.data[[colorFlag]])
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

  metadata2 <- data[data[[colorFlag]] == 1, ]
  metadata3 <- data[data[[colorFlag]] == -1, ]
  metadata4 <- rbind(metadata2, metadata3)

  metadata5 <- metadata4[metadata4[[expressionDirection]] == "up", ]
  metadata6 <- metadata4[metadata4[[expressionDirection]] == "down", ]

  meta.cor1 <- metadata5 %>%
    dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) %>%
    dplyr::summarize(n1 = paste(length(.data[[timeVariable]])), .groups = 'drop')

  p3 <- p2 + ggplot2::geom_text(data = meta.cor1, ggplot2::aes(x = 4, y = 9, label = n1),
                                colour = color3, inherit.aes = FALSE, parse = TRUE
  )

  meta.cor2 <- metadata6 %>%
    dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) %>%
    dplyr::summarize(n2 = paste(length(.data[[timeVariable]])), .groups = 'drop')

  p4 <- p3 + ggplot2::geom_text(data = meta.cor2, ggplot2::aes(x = -4, y = 9, label = n2),
                                colour = color1, inherit.aes = FALSE, parse = FALSE
  )

  # Return the created plot
  return(p4 + ggplot2::theme(legend.position = "bottom"))
}



setMethod("show",
          signature(object = "ClearScatterplot"),
          function(object) {
            object <- createPlot(object)
            print(object@plot)
          }
)
