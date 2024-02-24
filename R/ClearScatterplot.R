# Defines the S4 class 'ClearScatterplot' to represent scattered plot with specific data and plot configurations
#' ClearScatterplot: An S4 class for scatterplot visualization
#'
#' This class represents a scatterplot with specific data and plot configurations.
#'
#' @slot data_filepath Character vector: Path to data file.
#' @slot data Data.frame: data for the scatterplot.
#' @slot plot ANY: The ggplot2 plot object (initially NULL).
#' @name ClearScatterplot
#' @docType class
#' @importFrom methods setClass
#' @export
setClass("ClearScatterplot",
         slots = c(
           data = "data.frame",  # Store the dataframe directly
           plot = "ANY"              # Plot object (initially NULL)
         ),
         prototype = list(
           data = data.frame(),  # Default empty dataframe
           plot = NULL               # Default NULL plot
         )
)

# setClass("ClearScatterplot",
#          representation(
#            data = "data.frame",
#            plot = "ANY"
#          )
# )

#' @title Initialize a ClearScatterplot Object
#' @description Constructor method for the `ClearScatterplot` class.
#'   Initializes a `ClearScatterplot` object with specified data and plot configurations.
#' @param .Object The `ClearScatterplot` object to be initialized.
#' @param data_filepath Character string specifying the path to the data file.
#'   The data should contain information necessary for creating scatter plots.
#' @param timePointLevels Character vector specifying the levels for time points.
#'   Defaults to c("TP-1", "TP-2", "TP-3", "TP-4").
#' @param significanceColumn Character string specifying the name of the column in the data
#'   that contains the significant value. Example is "p", "q", "std". Default is "p"
#' @param expressionColumnName Character string specifying the name of the column in the data
#'   that contains expression values. Default is "log2fc".
#' @param highLog2fc Numeric value specifying the threshold for high log2 fold change.
#'   Default is 0.585.
#' @param lowLog2fc Numeric value specifying the threshold for low log2 fold change.
#'   Default is -0.585.
#' @param negLog10pValue Numeric value specifying the threshold for negative log10 p-values.
#'   Default is 1.301.
#' @return An initialized `ClearScatterplot` object with data loaded and ready for plot creation.
#' @export
#' @examples
#' \dontrun{
#'   data_path <- system.file("extdata", "example_data.csv", package = "YourPackage")
#'   csPlot <- new("ClearScatterplot", data_filepath = data_path)
#' }
# setMethod("initialize",
#           signature(.Object = "ClearScatterplot"),
#           function(.Object, data_filepath, timePointLevels = c("TP-1", "TP-2", "TP-3", "TP-4"),
#                    significanceColumn = "p", expressionColumnName = "log2fc",
#                    highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301) {
#             # Initialize a 'ClearScatterplot' object with the given data and default or specified parameters
#             .Object@data_filepath <- data_filepath # Path to data file
#             .Object@data <- process_data(data_filepath, significanceColumn = significanceColumn,
#                                                  expressionColumnName = expressionColumnName,
#                                                  highLog2fc = highLog2fc, lowLog2fc = lowLog2fc, negLog10pValue =  negLog10pValue)
#             .Object@plot <- NULL # Placeholder for the plot, to be created later
#             return(.Object) # Return the initialized object
#           }
# )
setMethod("initialize",
          signature(.Object = "ClearScatterplot"),
          function(.Object, data, timePointLevels = c("TP-1", "TP-2", "TP-3", "TP-4"),
                   significanceColumn = "p", expressionColumnName = "log2fc",
                   neglog10p = "neglog10p", color_flag = "color_flag",
                   highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301) {
            # Initialize a 'ClearScatterplot' object with the provided data and default or specified parameters
            .Object@data <- data # Directly assign the provided dataframe
            .Object@data <- process_data(data, significanceColumn = significanceColumn,
                                         expressionColumnName = expressionColumnName,
                                         neglog10p = neglog10p, color_flag = color_flag,
                                         highLog2fc = highLog2fc, lowLog2fc = lowLog2fc,
                                         negLog10pValue =  negLog10pValue)

            # validate_data(.Object@data, significanceColumn, expressionColumnName,
                              # highLog2fc, lowLog2fc, negLog10pValue) # Optional: Validate data structure
            .Object@plot <- NULL # Placeholder for the plot, to be created later
            return(.Object) # Return the initialized object
          }
)
# Define a generic method 'createPlot' to create plots for objects
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' @title Create and Assign a Plot to a ClearScatterplot Object
#' @description Generates a scatter plot based on the stored data and aesthetic parameters provided. This plot is then assigned to the `plot` slot of the `ClearScatterplot` object.
#' @param object A `ClearScatterplot` object for which the plot will be created.
#' @param color1 Character string specifying the color for high significance points. Default is "cornflowerblue".
#' @param color2 Character string specifying the color for points that do not meet the high or low significance thresholds. Default is "grey".
#' @param color3 Character string specifying the color for low significance points. Default is "indianred".
#' @param highLog2fc Numeric threshold for considering an expression value significantly high. Default is 0.585.
#' @param lowLog2fc Numeric threshold for considering an expression value significantly low. Default is -0.585.
#' @param expressionDirection Character string indicating the direction of regulation to be highlighted. Possible values are "up", "down", or "regulation" for both. Default is "regulation".
#' @param timeVariable Character string specifying the column in the data that represents the time variable. This parameter can be used to facet the plot by time points. Default is "reg_time_org".
#' @param xAxis Character string specifying the column name to be used as the x-axis variable. Default is "organ".
#' @param yAxis Character string specifying the column name to be used as the y-axis variable. Default is "timePoint".
#' @param changeMagnitude Character string specifying the column name to be used for the magnitude of change. This parameter allows dynamic specification of the x-axis variable in the plot. Default is "log2fc".
#' @param timeVariable Character string specifying the column in the data that represents the time variable. This parameter can be used to facet the plot by time points. Default is "reg_time_org".
#' @return The `ClearScatterplot` object with its `plot` slot updated to include the newly created plot.
#' @export
#' @examples
#' \dontrun{
#'   data_path <- system.file("extdata", "example_data.csv", package = "YourPackage")
#'   csPlot <- new("ClearScatterplot", data_filepath = data_path)
#'   csPlot <- createPlot(csPlot)
#'   # Now, csPlot contains a ggplot object in its plot slot that can be displayed using show(csPlot)
#' }
setMethod("createPlot",
          signature(object = "ClearScatterplot"),
          function(object, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                   highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                   timeVariable="reg_time_org", xAxis = "organ", yAxis="timePoint", changeMagnitude="log2fc") {
            # Create a plot based on the data and specified aesthetic parameters
            if (is.null(object@plot)) {
              View(object@data)
              object@plot <- create_plot(object@data, color1 = color1, color2 = color2, color3=color3,
                                         highLog2fc = highLog2fc, lowLog2fc = lowLog2fc,
                                         expressionDirection = expressionDirection,
                                         timeVariable = timeVariable, xAxis = xAxis, yAxis = yAxis,
                                         changeMagnitude=changeMagnitude)
            }
            return(object) # Return the object with its plot
          }
)

# Function to create the actual plot based on processed data and aesthetic choices
create_plot <- function(data, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                        highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                        timeVariable="reg_time_org", xAxis = "organ", yAxis="timePoint", changeMagnitude="log2fc") {
  colorFlag = "color_flag"
  facetFormula = paste0(xAxis, "~", yAxis)

  gp_obj <- ggplot2::ggplot(data = data, ggplot2::aes(
    x = data[[changeMagnitude]], y = neglog10p, color = as.factor(color_flag)
  )) +
    ggplot2::geom_point(alpha = 0.5, size = 1.75) +
    ggplot2::theme(legend.position = "none") + ggplot2::geom_jitter() +
    ggplot2::scale_color_manual(values = c("cornflowerblue", "grey", "indianred")) +
    ggplot2::labs(
      x = expression("log2 (fold change)"), y = expression("-log10 (p-value)")
    ) +
    #ggplot2::facet_grid( ~ organ+timePoint, space = "free") +
    ggplot2::facet_grid(facetFormula, space = "free") +
    ggplot2::theme_bw() + ggplot2::theme(
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

  data2 <- data[data[[colorFlag]] == 1, ]
  data3 <- data[data[[colorFlag]] == -1, ]
  data4 <- rbind(data2, data3)

  data5 <- data4[data4[[expressionDirection]] == "up", ]
  data6 <- data4[data4[[expressionDirection]] == "down", ]

  meta.cor1 <- dplyr::summarize(
    dplyr::group_by(data5, data5[[xAxis]], data5[[yAxis]]),
    n1 = paste(length(.data[[timeVariable]])),
    .groups = 'drop'
  )

  p3 <- p2 + ggplot2::geom_text(data = meta.cor1, ggplot2::aes(x = 4, y = 9, label = n1),
                                colour = color3, inherit.aes = FALSE, parse = TRUE
  )

  meta.cor2 <- dplyr::summarize(
    dplyr::group_by(data6, data6[[xAxis]], data6[[yAxis]]),
    n2 = paste(length(.data[[timeVariable]])),
    .groups = 'drop'
  )


  p4 <- p3 + ggplot2::geom_text(data = meta.cor2, ggplot2::aes(x = -4, y = 9, label = n2),
                                colour = color1, inherit.aes = FALSE, parse = FALSE
  )

  # Return the created plot
  return(p4 + ggplot2::theme(legend.position = "bottom"))
}


#' @title Display the Scatterplot Stored in a ClearScatterplot Object
#' @description This method invokes the plotting functionality to generate and display the scatterplot based on the data and configuration stored within the `ClearScatterplot` object. If the plot has not been created yet, it first calls the `createPlot` method to generate the plot and then displays it.
#' @param object A `ClearScatterplot` object containing the necessary data and configuration for generating the scatterplot.
#' @return Invisible `ClearScatterplot` object. The primary purpose of this method is the side effect of displaying the plot, not returning a value. However, the object itself is returned invisibly for potential further manipulation or inspection.
#' @export
#' @examples
#' \dontrun{
#'   data_path <- system.file("extdata", "example_data.csv", package = "YourPackage")
#'   csPlot <- new("ClearScatterplot", data_filepath = data_path)
#'   # Assuming createPlot method has been properly defined and exported
#'   show(csPlot)
#'   # This will display the scatterplot based on the data
#' }
setMethod("show",
          signature(object = "ClearScatterplot"),
          function(object) {
            object <- createPlot(object)
            print(object@plot)
          }
)

# Function to process data for 'ClearScatterplot' objects
process_data <- function(data, timePointLevels, significanceColumn = "p",
                         expressionColumnName = "log2fc", neglog10p = "neglog10p",
                         color_flag = "color_flag",
                         highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301) {
  # data <- read.csv(filepath)
  #. data$neglog10p = -log10(.data[[significanceColumn]])
  #View(data)

  .data$color_flag <-
    ifelse(
      data$log2fc > highLog2fc &
        data$eglog10p > negLog10pValue, 1, ifelse(data$log2fc < lowLog2fc &
                                                         data$neglog10p > negLog10pValue,-1, 0)
    )

  .data$timePoint <- factor(
    data$timePoint, levels = timePointLevels
  )
  #View(data)
  return(data)
}

