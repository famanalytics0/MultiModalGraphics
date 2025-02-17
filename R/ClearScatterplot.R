# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2","neglog10p"))

#' ClearScatterplot: A Class for Enhanced Scatterplot
#'
#' A class for creating and managing scatterplot objects with enhanced visualization features.
#' @slot data A data frame containing the data to be plotted.
#' @slot plot An object storing the ggplot representation of the data.
#' @name ClearScatterplot
#' @docType class
#' @importFrom MultiAssayExperiment experiments colData assay
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix
#' @importFrom ggplot2 aes geom_point geom_text geom_jitter facet_grid
#'   theme scale_fill_hue scale_color_manual labs element_text element_rect
#' @export ClearScatterplot
#' @exportMethod ClearScatterplot
#' @exportMethod createPlot
#' @exportMethod show
setClass("ClearScatterplot",
         slots = c(
           data = "data.frame", # Slot for storing the data frame
           plot = "ANY" # Slot for storing the ggplot object
         )
)

setGeneric("ClearScatterplot", function(data, ...) {
  standardGeneric("ClearScatterplot")
})

#' ClearScatterplot Constructor
#'
#' Creates an instance of the ClearScatterplot class.
#' @param data A data frame containing the plot data.
#' @param logFoldChange The name of the column containing expression values.
#' @param negativeLogPValue The column containing the negative log p-values.
#' @param highLog2fc Threshold for high log2 fold change values.
#' @param lowLog2fc Threshold for low log2 fold change values.
#' @param negLog10pValue Threshold for -log10 p-value.
#' @param timePointColumn The name of the column containing time point info.
#' @param timePointLevels The levels for the time point column, if any.
#' @return An object of class ClearScatterplot.
#' @export
#' @examples
#' plotdata <- get_clear_scatterplot_df()
#' scatterplotObject <- ClearScatterplot(
#'   data = plotdata,
#'   logFoldChange = "log2fc",
#'   timePointColumn = "timePoint",
#'   timePointLevels = c("T10R1", "T5R1")
#' )
ClearScatterplot <- function(data,
                             logFoldChange = "log2fc",
                             negativeLogPValue = "negLog10p",
                             highLog2fc = 0.585,
                             lowLog2fc = -0.585,
                             negLog10pValue = 1.301,
                             timePointColumn = "timePoint",
                             timePointLevels = NULL) {
  if (!is.data.frame(data)) {
    stop("Data must be a data frame.")
  }

  # Data preprocessing
  data$color_flag <- ifelse(
    data[[logFoldChange]] > highLog2fc &
      data[[negativeLogPValue]] > negLog10pValue,
    1,
    ifelse(data[[logFoldChange]] < lowLog2fc &
             data[[negativeLogPValue]] > negLog10pValue, -1, 0)
  )

  # Factorize the timePoint column based on specified levels if provided
  if (!is.null(timePointLevels) &&
      !is.null(data[[timePointColumn]])) {
    data[[timePointColumn]] <-
      factor(data[[timePointColumn]], levels = timePointLevels)
  }

  # Create the object
  obj <- new("ClearScatterplot", data = data)

  return(obj)
}



#' Create a ClearScatterplot Object from MultiAssayExperiment
#'
#' Constructor function for creating an instance of the ClearScatterplot class from a MultiAssayExperiment object.
#' This function performs differential expression analysis using limma.
#' @param mae A MultiAssayExperiment object containing the assay data.
#' @param assayName The name of the assay to use for plotting.
#' @param logFoldChange The name of the column containing expression values.
#' @param negativeLogPValue The name of the column containing the negative log p-values.
#' @param highLog2fc Threshold for high log2 fold change values.
#' @param lowLog2fc Threshold for low log2 fold change values.
#' @param negLog10pValue Threshold for -log10 p-value.
#' @param timepoint The name of the column containing time point information.
#' @param sampleType The name of the column containing sample type information (e.g., tissue or organ).
#' @return An object of class ClearScatterplot.
#' @examples
#' # Parameters
#' num_genes <- 1000
#' num_samples <- 1000  # Total number of samples
#' conditions <- list(
#'   Condition1 = c("Cond1A", "Cond1B"),
#'   Condition2 = c("Cond2A", "Cond2B"),
#'   Condition3 = c("Cond3A", "Cond3B")
#' )
#'
#' # Function to generate expression data
#' generate_expression_data <- function(num_genes, num_samples) {
#'   matrix(rnorm(num_genes * num_samples, mean = 6, sd = 2), nrow = num_genes, ncol = num_samples)
#' }
#'
#' # Initialize lists
#' expression_data_list <- list()
#' metadata_list <- list()
#'
#' # Generate samples for each combination of conditions
#' all_combinations <- expand.grid(conditions)
#' num_combinations <- nrow(all_combinations)
#' samples_per_combination <- num_samples %/% num_combinations
#' remaining_samples <- num_samples %% num_combinations
#'
#' for (i in 1:num_combinations) {
#'   combination <- all_combinations[i, ]
#'   combination_name <- paste(combination, collapse = "_")
#'
#'   # Calculate the number of samples for this combination
#'   num_samples_for_this_combination <- samples_per_combination + ifelse(i <= remaining_samples, 1, 0)
#'
#'   expression_data <- generate_expression_data(num_genes, num_samples_for_this_combination)
#'   colnames(expression_data) <- paste0(combination_name, "_Sample", 1:num_samples_for_this_combination)
#'   rownames(expression_data) <- paste0("Gene", 1:num_genes)
#'
#'   sample_metadata <- data.frame(
#'     SampleID = colnames(expression_data),
#'     Group = rep(c("Control", "Treatment"), length.out = num_samples_for_this_combination),
#'     stringsAsFactors = FALSE,
#'     TimePoint = rep(c("T1", "T1", "T2", "T2"), length.out = num_samples_for_this_combination),  # Example time points
#'     SampleType = rep(c("Tissue", "Organ"), length.out = num_samples_for_this_combination)  # Example sample types
#'   )
#'
#'   for (cond in names(combination)) {
#'     sample_metadata[[cond]] <- combination[[cond]]
#'   }
#'
#'   expression_data_list[[combination_name]] <- expression_data
#'   metadata_list[[combination_name]] <- sample_metadata
#' }
#'
#' # Combine metadata into a single data frame
#' combined_metadata <- do.call(rbind, metadata_list)
#' rownames(combined_metadata) <- combined_metadata$SampleID
#'
#' # Select only the relevant columns for the final metadata
#' final_metadata <- combined_metadata[, c("SampleID", "TimePoint", "SampleType", "Group", "Condition1", "Condition2", "Condition3")]
#'
#' # Create the MultiAssayExperiment object
#' mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = expression_data_list, colData = final_metadata)
#'
#' # Check available assay names
#' assayNames <- names(MultiAssayExperiment::experiments(mae))
#' assayName <- assayNames[1]
#'
#' # Assuming ClearScatterplot_MAE is already defined in your script
#' scatterplotObject <- ClearScatterplot_MAE(
#'   mae = mae,
#'   assayName = assayName,
#'   timepoint = "TimePoint",
#'   sampleType = "SampleType"
#' )
#'
#' # Assuming createPlot is already defined in your script
#' # Generate and display the plot
#' scatterplotObject <- createPlot(
#'   scatterplotObject,
#'   color1 = "cornflowerblue",
#'   color2 = "grey",
#'   color3 = "indianred",
#'   highLog2fc = 0.585,
#'   lowLog2fc = -0.585,
#'   negLog10pValue = 1.301
#' )
#'
#' # Show the plot
#' print(scatterplotObject@plot)
#' @export
ClearScatterplot_MAE <- function(mae, assayName = NULL, logFoldChange = "log2fc",
                                 negativeLogPValue = "negLog10p", highLog2fc = 0.585,
                                 lowLog2fc = -0.585, negLog10pValue = 1.301,
                                 timepoint = "TimePoint", sampleType = "SampleType") {
    if (!inherits(mae, "MultiAssayExperiment")) {
    stop("Data must be a MultiAssayExperiment object.")
  }

  if (is.null(assayName) || !assayName %in% names(MultiAssayExperiment::experiments(mae))) {
    stop("Please specify a valid assayName from the MultiAssayExperiment object.")
  }

  se <- MultiAssayExperiment::experiments(mae)[[assayName]]
  col_data <- MultiAssayExperiment::colData(mae)
  expr_data <- MultiAssayExperiment::assay(se)

  # Ensure columns match and are in the correct order
  sample_names <- colnames(expr_data)
  col_data <- col_data[rownames(col_data) %in% sample_names, ]
  col_data <- col_data[match(sample_names, rownames(col_data)), ]

  # Differential expression analysis using limma
  design <- model.matrix(~ Group, data = col_data)
  fit <- limma::lmFit(expr_data, design)
  fit <- limma::eBayes(fit)
  topTable <- limma::topTable(fit, adjust.method = "BH", number = Inf)

  # Format the results
  data <- data.frame(
    log2fc = topTable$logFC,
    negLog10p = -log10(topTable$P.Value),
    regulation = ifelse (topTable$logFC > 0, "up", "down"),
    organ = colData(mae)[[sampleType]],
    timePoint = colData(mae)[[timepoint]],
    reg_time_org = ifelse (topTable$logFC > 0, "up", "down"),
    row.names = rownames(topTable)
  )

  ClearScatterplot(data, logFoldChange, negativeLogPValue, highLog2fc, lowLog2fc, negLog10pValue)
}





# Define a generic method 'createPlot' to create plots for objects
#' Define a Generic Method 'createPlot'
#'
#' This generic function is designed to create plots for objects of various classes.
#' Specific methods should be defined for different classes to generate the appropriate plots.
#'
#' @param object The object for which the plot is to be created. The class of this object
#'   will determine which specific method is used.
#' @param ... Additional arguments to be passed to the specific method for creating the plot.
#'
#' @return The result of the specific method for creating the plot, typically a plot object.
#' @rdname createPlot
#' @export
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
#' @param negLog10pValue Threshold for -log10 p-value.
#' @param negativeLogPValue The name of the column containing the negative log p-values.
#' @param timeVariable The variable representing time.
#' @param xAxis The x-axis values.
#' @param yAxis The y-axis values.
#' @return The ClearScatterplot object with the plot updated.
#' @rdname createPlot
#' @export
#' @examples
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
setMethod(
  "createPlot",
  signature(object = "ClearScatterplot"),
  function(object,
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
           yAxis = "timePoint") {
    # Create a plot based on the data and specified aesthetic parameters
    if (is.null(object@plot)) {
      object@plot <-
        create_plot(
          object@data,
          color1 = color1,
          color2 = color2,
          color3 = color3,
          highLog2fc = highLog2fc,
          lowLog2fc = lowLog2fc,
          expressionDirection = expressionDirection,
          negativeLogPValue = negativeLogPValue,
          timeVariable = timeVariable,
          xAxis = xAxis,
          yAxis = yAxis
        )
    }
    return(object) # Return the object with its plot
  }
)

#' Create the Actual Plot
#'
#' This function creates the actual plot based on processed data and aesthetic choices.
#' @param data The data frame containing the plot data.
#' @param color1 Color for one category of data points.
#' @param color2 Color for another category of data points.
#' @param color3 Color for a third category of data points.
#' @param highLog2fc Threshold for high log2 fold change values.
#' @param lowLog2fc Threshold for low log2 fold change values.
#' @param expressionDirection Direction of gene expression.
#' @param negativeLogPValue The name of the column containing the negative log p-values.
#' @param logFoldChange The name of the column containing log fold change values.
#' @param timeVariable The variable representing time.
#' @param xAxis The name of the x-axis variable.
#' @param yAxis The name of the y-axis variable.
#' @return The ggplot object.
#' @export
#' @examples
#' plotdata <- get_clear_scatterplot_df()
#' p <- create_plot(
#'   plotdata,
#'   color1 = "cornflowerblue",
#'   color2 = "grey",
#'   color3 = "indianred",
#'   highLog2fc = 0.585,
#'   lowLog2fc = -0.585,
#'   expressionDirection = "regulation",
#'   negativeLogPValue = "negLog10p",
#'   logFoldChange = "log2fc",
#'   timeVariable = "reg_time_org",
#'   xAxis = "organ",
#'   yAxis = "timePoint"
#' )
create_plot <- function(data,
                        color1 = "cornflowerblue",
                        color2 = "grey",
                        color3 = "indianred",
                        highLog2fc = 0.585,
                        lowLog2fc = -0.585,
                        expressionDirection = "regulation",
                        negativeLogPValue = "negLog10p",
                        logFoldChange = "log2fc",
                        timeVariable = "reg_time_org",
                        xAxis = "organ",
                        yAxis = "timePoint") {
  colorFlag <- "color_flag"
  facetFormula <- paste0(xAxis, "~", yAxis)

  gp_obj <- ggplot(data = data, ggplot2::aes(
    x = .data[[logFoldChange]],
    y = .data[[negativeLogPValue]],
    color = as.factor(.data[[colorFlag]])
  )) +
    ggplot2::geom_point(alpha = 0.5, size = 1.75) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_jitter() +
    ggplot2::scale_color_manual(values = c(color1, color2, color3)) +
    ggplot2::labs(
      x = expression("log2 (fold change)"), y = expression("-log10 (p-value)")
    ) +
    ggplot2::facet_grid(facetFormula, space = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12)
    )

  gp_obj2 <- gp_obj + ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(
      size = 12,
      face = "bold"
    ),
    title = ggplot2::element_text(size = 12, face = "bold"),
    strip.text.x = ggplot2::element_text(
      size = 12,
      face = "bold",
      colour = "black",
      angle = 0
    )
  )
  gp_obj3 <- gp_obj2 + ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),
    axis.title.y = ggplot2::element_text(
      size = 12,
      face = "bold"
    ),
    title = ggplot2::element_text(size = 12, face = "bold"),
    strip.text.y = ggplot2::element_text(
      size = 12,
      face = "bold",
      colour = "black",
      angle = 90
    )
  )
  p <- gp_obj3 + ggplot2::theme(legend.position = "right")

  p1 <-
    p + ggplot2::theme(strip.background = ggplot2::element_rect(
      fill = "white", color = color2, linewidth = 1
    ))

  p2 <-
    p1 + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::scale_color_manual(
      name = "log2 (fold change)",
      values = c(color1, color2, color3),
      labels = c(
        paste("<", lowLog2fc),
        paste(lowLog2fc, "to", highLog2fc),
        paste(">", highLog2fc)
      )
    ) + ggplot2::scale_fill_hue()

  metadata2 <- data[data[[colorFlag]] == 1, ]
  metadata3 <- data[data[[colorFlag]] == -1, ]
  metadata4 <- rbind(metadata2, metadata3)

  metadata5 <- metadata4[metadata4[[expressionDirection]] == "up", ]
  metadata6 <- metadata4[metadata4[[expressionDirection]] == "down", ]

  meta.cor1 <- metadata5 %>%
    dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) %>%
    dplyr::summarize(n1 = paste(length(.data[[timeVariable]])), .groups = "drop")

  p3 <-
    p2 + ggplot2::geom_text(
      data = meta.cor1,
      ggplot2::aes(x = 4, y = 9, label = n1),
      colour = color3,
      inherit.aes = FALSE,
      parse = TRUE
    )

  meta.cor2 <- metadata6 %>%
    dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) %>%
    dplyr::summarize(n2 = paste(length(.data[[timeVariable]])), .groups = "drop")

  p4 <-
    p3 + ggplot2::geom_text(
      data = meta.cor2,
      ggplot2::aes(x = -4, y = 9, label = n2),
      colour = color1,
      inherit.aes = FALSE,
      parse = FALSE
    )

  # Return the created plot
  return(p4 + ggplot2::theme(legend.position = "bottom"))
}

#' Display Method for ClearScatterplot
#'
#' This method displays the scatterplot stored in the ClearScatterplot object.
#' It invokes the `print` method on the `plot` object, which must be a ggplot object.
#' @param object A ClearScatterplot object.
#' @return Prints the ggplot object stored in the plot slot of the ClearScatterplot.
#' The function itself returns invisible `NULL` to avoid additional console output.
#' @export
#' @examples
#' # Create example data
#' data <- data.frame(log2fc = rnorm(100), negLog10p = -log10(runif(100)))
#' scatterplotObject <- ClearScatterplot(data)
#' 
#' # Display the plot
#' show(scatterplotObject)
setMethod(
  "show",
  signature(object = "ClearScatterplot"),
  function(object) {
    if (is.null(object@plot)) {
      object <- createPlot(object) # Update the plot only if it's null
    }
    # Display the plot
    print(object@plot)
    # Optionally return the object invisibly for further chaining if needed
    invisible(object)
  }
)
