# Suppress warnings for global variables
utils::globalVariables(c("assays", "setNames", "iClusterPlus"))

#' InformativeHeatmap: A Class for Enhanced Heatmaps
#'
#' This class encapsulates the ComplexHeatmap package's Heatmap object,
#' allowing for extended functionality and customizations.
#'
#' @slot heatmap Heatmap object from ComplexHeatmap package.
#' @slot params List of parameters used to construct the heatmap.
#' @name InformativeHeatmap
#' @docType class
#' @importFrom ComplexHeatmap Heatmap
#' @export InformativeHeatmap
#' @exportMethod InformativeHeatmap
#' @exportMethod updateLayerFun
#' @exportMethod getHeatmapObject

setClass("InformativeHeatmap",
         slots = c(
           heatmap = "Heatmap",  # The actual heatmap object
           params = "list"  # Slot to store the constructor parameters
         ))

setGeneric("InformativeHeatmap", function(data, ...) {
  standardGeneric("InformativeHeatmap")
})

#' @rdname InformativeHeatmap
#' @description
#' Constructor for the `InformativeHeatmap` class. This function initializes
#' an `InformativeHeatmap` object with the given data and additional parameters
#' for heatmap customization, including visual customization based on significance levels.
#'
#' @param data A numeric matrix or data frame to be visualized as a heatmap.
#' @param pch_val Plotting character (pch) value for points in the heatmap.
#'   Default is 16.
#' @param unit_val Size of the points in the heatmap. Specified in 'mm'.
#'   Default is 1.
#' @param significant_color Color to be used for points representing significant values.
#'   Default is "black".
#' @param trending_color Color to be used for points representing trending values
#'   (significant but less so than those meeting the `significant_pvalue` criterion).
#'   Default is "yellow".
#' @param significant_pvalue P-value threshold for significance.
#'   Points with p-values below this threshold will be colored with `significant_color`.
#'   Default is 0.05.
#' @param trending_pvalue P-value threshold for trending significance.
#'   Points with p-values above `significant_pvalue` and below this threshold
#'   will be colored with `trending_color`. Default is 0.1.
#' @param ... Additional parameters passed to `ComplexHeatmap::Heatmap`.
#' @return An object of class `InformativeHeatmap` with the heatmap initialized
#'   and customized according to the provided parameters.
#' @examples
#' \donttest{
#'   data <- matrix(rnorm(100), ncol = 10)
#'   heatmap <- InformativeHeatmap(data, pch_val = 20, unit_val = 2,
#'                                 significant_color = "red",
#'                                 trending_color = "blue",
#'                                 significant_pvalue = 0.05,
#'                                 trending_pvalue = 0.1)
#' }
#' @export

setMethod("InformativeHeatmap", "ANY", function(data,
                                                pch_val = 16,
                                                unit_val = 4,
                                                significant_color = "black",
                                                trending_color = "yellow",
                                                significant_pvalue = 0.05,
                                                trending_pvalue = 0.1,
                                                ...) {
  # First, check if 'data' is a matrix
  if (!is.matrix(data)) {
    stop("Data must be a matrix.")
  }

  # Second, check if the matrix contains only numeric values
  if (!all(is.numeric(data))) {
    stop("Data matrix must contain only numeric values.")
  }

  if (!requireNamespace("ComplexHeatmap", quietly <- TRUE)) {
    stop(
      "ComplexHeatmap is required for creating an InformativeHeatmap object. ",
      "Please install it using BiocManager::install('ComplexHeatmap')."
    )
  }
  # Capture all additional arguments passed to the function
  params_list <- list(...)

  # Initialize var for potential custom layer function based on small_map_pv
  custom_layer_fun <- NULL

  # Now include layer_fun directly in the params, if not already provided
  if (!is.null(params_list$significance_level)) {
    # Store the small_map_pv data
    significance_level <- params_list$significance_level

    custom_layer_fun <- function(j, i, x, y, w, h, fill) {
      ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
      for (ir in seq_len(nrow(ind_mat))) {
        for (ic in seq_len(ncol(ind_mat))) {
          ind <- ind_mat[ir, ic]  # previous column
          v <- significance_level[i[ind], j[ind]]
          grid::grid.points(x[ind],
                      y[ind],
                      pch <-
                        pch_val,
                      gp = grid::gpar(col = ifelse(
                        v < significant_pvalue,
                        significant_color,
                        ifelse(
                          v >= significant_pvalue && v < trending_pvalue,
                          trending_color,
                          NA
                        )
                      )),
                      size <- grid::unit(unit_val, "mm"))
        }
      }
    }
    # Remove small_map_pv from params_list to prevent the unused argument error
    params_list$significance_level <- NULL
  }

  # If custom_layer_fun is defined, update layer_fun in params_list
  if (!is.null(custom_layer_fun)) {
    params_list$layer_fun <- custom_layer_fun
  }

  # Create the heatmap object with the remaining and possibly modified parameters
  heatmap_obj <- do.call(ComplexHeatmap::Heatmap, c(list(data), params_list))

  # Construct and return the InformativeHeatmap object
  methods::new("InformativeHeatmap", heatmap = heatmap_obj, params = params_list)

})


#' Define a Generic Method 'updateLayerFun'
#'
#' This generic function is for updating the layer function of an existing
#' InformativeHeatmap objects.
#'
#' @param x An `InformativeHeatmap` object whose layer function is to be updated.
#' @param layer_fun A function that defines the new layer to be applied to the heatmap.
#'   The class of this object will determine which specific method is used.
#'
#' @return Returns an updated `InformativeHeatmap` object with the new layer function applied.
#' @export
setGeneric("updateLayerFun", function(x, layer_fun) {
  standardGeneric("updateLayerFun")
})

#' Update Layer Function of an InformativeHeatmap
#'
#' This method allows updating the layer function of an existing `InformativeHeatmap` object.
#' It requires the `ComplexHeatmap` package to recreate the heatmap with the new layer function.
#' If `ComplexHeatmap` is not installed, it will stop and prompt the user to install `ComplexHeatmap`.
#'
#' @param x An `InformativeHeatmap` object whose layer function is to be updated.
#' @param layer_fun A function that defines the new layer to be applied to the heatmap.
#'   This function should be compatible with the layering system of `ComplexHeatmap`.
#' @return Returns an updated `InformativeHeatmap` object with the new layer function applied.
#' @examples
#' \donttest{
#'   # Assume `ih` is an existing InformativeHeatmap object
#'   # Define a new layer function
#'   new_layer_fun <- function(j, i, x, y, w, h, fill) {
#'     grid::grid.points(x, y, pch = 16, size = unit(2, "mm"), gp = grid::gpar(col = "red"))
#'   }
#'
#'   # Update the layer function of the heatmap
#'   ih <- updateLayerFun(ih, new_layer_fun)
#' }
#' @export
setMethod("updateLayerFun", "InformativeHeatmap", function(x, layer_fun) {
  if (!requireNamespace("ComplexHeatmap", quietly <- TRUE)) {
    stop(
      "ComplexHeatmap is required to update the layer function in an InformativeHeatmap object. ",
      "Please install it using BiocManager::install('ComplexHeatmap')."
    )
  }
  # Retrieve the stored parameters
  params <- x@params

  # Update or add the layer_fun parameter
  params$layer_fun <- layer_fun

  # Ensure the matrix data is passed as the first argument to Heatmap
  # and the rest of the parameters are correctly structured for do.call
  args <- c(list(x@heatmap@matrix), params)

  # Recreate the heatmap with the updated parameters
  new_heatmap <- do.call(ComplexHeatmap::Heatmap, args)

  # Update the InformativeHeatmap object with the new heatmap
  x@heatmap <- new_heatmap

  # Also, update the stored parameters in case of further modifications
  x@params <- params

  return(x)
})

# setGeneric("InformativeHeatmapFromMAE", function(mae, ...) {
#   standardGeneric("InformativeHeatmapFromMAE")
# })


#' Create an InformativeHeatmap from a MultiAssayExperiment Object
#'
#' This method takes a `MultiAssayExperiment` object, clusters the samples using `iClusterPlus`,
#' compares the clusters using `limma`, and creates an `InformativeHeatmap` with a custom `layer_fun`
#' based on significance levels from the cluster comparison.
#'
#' @param mae A `MultiAssayExperiment` object containing multiple assays.
#' @param significant_pvalue P-value threshold for significance.
#'   Points with p-values below this threshold will be colored with `significant_color`.
#'   Default is 0.05.
#' @param trending_pvalue P-value threshold for trending significance.
#'   Points with p-values above `significant_pvalue` and below this threshold
#'   will be colored with `trending_color`. Default is 0.1.
#' @param significant_color Color to be used for points representing significant values.
#'   Default is "black".
#' @param trending_color Color to be used for points representing trending values.
#'   Default is "yellow".
#' @param pch_val Plotting character (pch) value for points in the heatmap. Default is 16.
#' @param unit_val Size of the points in the heatmap. Specified in 'mm'. Default is 1.
#' @param K Number of clusters for `iClusterPlus`. Default is 3.
#' @param lambda Regularization parameter for `iClusterPlus`. Default is 0.2.
#' @param coef Compare cluster group to be used as a contrast for `limma` analysis.
#' @param ... Additional parameters passed to `ComplexHeatmap::Heatmap`.
#' @return An object of class `InformativeHeatmap`.
#' @export
InformativeHeatmapFromMAE <- function(mae,
                                      significant_pvalue = 0.05,
                                      trending_pvalue = 0.1,
                                      significant_color = "black",
                                      trending_color = "yellow",
                                      pch_val = 16,
                                      unit_val = 4,
                                      K = 3,  # Number of clusters for iClusterPlus
                                      lambda = 0.2,  # Regularization for iClusterPlus
                                      coef = 2,  # Contrast to compare cluster groups in limma
                                      ...) {

  # Convert MultiAssayExperiment to a list of matrices for iClusterPlus
  data_list <- lapply(assays(mae), function(assay) {
    assay_matrix <- as.matrix(assay)
    mode(assay_matrix) <- "numeric"
    return(t(assay_matrix))  # Transpose so that samples are rows and features are columns
  })

  # Run iClusterPlus clustering
  set.seed(123)
  dt_list <- stats::setNames(data_list, paste0("dt", seq_along(data_list)))

  fit <- do.call(iClusterPlus, c(dt_list, list(
    type = rep("gaussian", length(data_list)),
    K = K,  # Number of clusters
    lambda = rep(lambda, length(data_list)),
    maxiter = 20
  )))

  # Get the clustering results
  clusters <- factor(fit$clusters)

  # Create design matrix for limma analysis based on clusters
  design <- model.matrix(~ 0 + clusters)
  colnames(design) <- levels(clusters)

  # Perform limma analysis on each assay to compare the cluster groups
  limma_results <- lapply(names(assays(mae)), function(assay_name) {
    assay_data <- assays(mae)[[assay_name]]
    fit_limma <- lmFit(assay_data, design)
    fit_limma <- eBayes(fit_limma)
    topTable_limma <- topTable(fit_limma, coef = coef, number = Inf)

    # Extract logFC and p-values
    logFC <- topTable_limma$logFC
    p_values <- topTable_limma$P.Value

    return(list(logFC = logFC, p_values = p_values))
  })
  names(limma_results) <- names(assays(mae))

  # Combine logFCs across assays for heatmap
  combined_logFC <- do.call(rbind, lapply(limma_results, function(result) {
    return(result$logFC)
  }))

  # Combine p-values across assays
  combined_pvalues <- do.call(rbind, lapply(limma_results, function(result) {
    return(result$p_values)
  }))

  # Ensure that combined_pvalues and combined_logFC have the same number of rows
  if (nrow(combined_logFC) != nrow(combined_pvalues)) {
    stop("Mismatch between the number of rows in logFC and p-values")
  }

  # Transpose combined data matrix for heatmap
  heatmap_data <- t(combined_logFC)

  # Define a custom layer_fun based on p-values
  custom_layer_fun <- function(j, i, x, y, w, h, fill) {
    ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
    n_rows <- nrow(combined_pvalues)
    n_cols <- ncol(combined_pvalues)

    for (ir in seq_len(nrow(ind_mat))) {
      for (ic in seq_len(ncol(ind_mat))) {
        ind <- ind_mat[ir, ic]
        # Make sure to use the right index to access combined_pvalues and combined_logFC
        if (ir <= n_rows && ic <= n_cols) {
          p_value <- combined_pvalues[ir, ic]

          # Determine color based on p-value
          col_value <- ifelse(
            p_value < significant_pvalue,
            significant_color,
            ifelse(
              p_value >= significant_pvalue && p_value < trending_pvalue,
              trending_color,
              NA
            )
          )

          # Only draw points for non-NA values
          if (!is.na(col_value)) {
            grid::grid.points(
              x[ind], y[ind],
              pch = pch_val,
              gp = grid::gpar(col = col_value),
              size = grid::unit(unit_val, "mm")
            )
          }
        } else {
          warning("Index out of bounds in combined_pvalues")
        }
      }
    }
  }

  # Create and return the InformativeHeatmap object with the custom layer_fun
  return(InformativeHeatmap(heatmap_data, layer_fun = custom_layer_fun, ...))
}

#' Define a Generic Method 'getHeatmapObject'
#'
#' This generic function is designed to retrieve the Heatmap object from
#' InformativeHeatmap objects. Specific methods should be defined for
#' different classes to extract the Heatmap object as needed.
#'
#' @param x The object from which the Heatmap object is to be retrieved.
#'   The class of this object will determine which specific method is used.
#'
#' @return The result of the specific method for retrieving the Heatmap object,
#'   typically a Heatmap object.
#' @export
setGeneric("getHeatmapObject", function(x) {
  standardGeneric("getHeatmapObject")
})

#' Retrieve the Heatmap Object from an InformativeHeatmap
#'
#' This method extracts the underlying Heatmap object stored within an
#' `InformativeHeatmap` object. It allows direct access to the `Heatmap` object
#'  for further manipulation or inspection using `ComplexHeatmap` package
#'  functionalities. Note that the `ComplexHeatmap` package is required
#' to fully utilize the returned Heatmap object.
#'
#' @param x An `InformativeHeatmap` object from which the Heatmap object is to
#' be retrieved.
#'
#' @return A Heatmap object from the `ComplexHeatmap` package.
#'
#' @examples
#' \donttest{
#'   # Assume `ih` is an existing InformativeHeatmap object
#'   heatmap_obj <- getHeatmapObject(ih)
#'   # Now `heatmap_obj` can be used directly with ComplexHeatmap functions
#' }
#'
#' @export
setMethod("getHeatmapObject", "InformativeHeatmap", function(x) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop(
      "ComplexHeatmap is required to retrieve the Heatmap object from an
      InformativeHeatmap object. ",
      "Please install it using BiocManager::install('ComplexHeatmap')."
    )
  }
  return(x@heatmap)
})





