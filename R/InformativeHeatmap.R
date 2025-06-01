# Suppress notes about global variables
utils::globalVariables(c("assays", "setNames", "iClusterPlus"))

#’ InformativeHeatmap: A Class for Enhanced Heatmaps
#’
#’ This class encapsulates a ComplexHeatmap::Heatmap object, allowing for extended
#’ functionality and customizations (e.g., overlaying points colored by p-value).
#’
#’ @slot heatmap Heatmap object from ComplexHeatmap package.
#’ @slot params List of parameters used to construct the heatmap.
#’ @name InformativeHeatmap
#’ @docType class
#’ @importFrom ComplexHeatmap Heatmap restore_matrix
#’ @importFrom grid grid.points gpar unit
#’ @importFrom stats model.matrix lm
#’ @exportClass InformativeHeatmap

setClass(
  "InformativeHeatmap",
  slots = c(
    heatmap = "Heatmap",
    params  = "list"
  ),
  validity = function(object) {
    if (!inherits(object@heatmap, "Heatmap")) {
      return("Slot ‘heatmap’ must be a ComplexHeatmap::Heatmap object.")
    }
    if (!is.list(object@params)) {
      return("Slot ‘params’ must be a list.")
    }
    TRUE
  }
)

#’ @rdname InformativeHeatmap
#’ @exportMethod InformativeHeatmap
#’ @export
setGeneric("InformativeHeatmap", function(data, ...) standardGeneric("InformativeHeatmap"))

#’ Constructor for InformativeHeatmap
#’
#’ Builds an InformativeHeatmap object from a numeric matrix, with optional overlay
#’ of points colored by significance/trending thresholds.
#’
#’ @param data A numeric matrix (rows = features, columns = samples).
#’ @param pch_val Integer or single‐length numeric. Plotting character for overlay points (default = 16).
#’ @param unit_val Numeric. Size of overlay points in mm (default = 4).
#’ @param significant_color Color for points with p < significant_pvalue (default = "black").
#’ @param trending_color Color for points with significant_pvalue ≤ p < trending_pvalue (default = "yellow").
#’ @param significant_pvalue Numeric ∈ (0,1). Threshold for “significant” (default = 0.05).
#’ @param trending_pvalue Numeric ∈ (0,1). Threshold for “trending” (default = 0.1).
#’ @param ... Additional arguments passed to ComplexHeatmap::Heatmap (e.g. name, col, cluster_rows, etc.).
#’
#’ @return An object of class InformativeHeatmap.
#’ @examples
#’ \dontrun{
#’   mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
#’   # Create a matrix of random p‐values (same dims as mat)
#’   pmat <- matrix(runif(200), nrow = 20, ncol = 10)
#’   ih <- InformativeHeatmap(
#’     mat,
#’     significance_level = pmat,
#’     pch_val            = 20,
#’     unit_val           = 2,
#’     significant_color  = "red",
#’     trending_color     = "blue",
#’     significant_pvalue = 0.01,
#’     trending_pvalue    = 0.05,
#’     name               = "Expression",
#’     cluster_rows       = TRUE,
#’     cluster_columns    = TRUE
#’   )
#’   hm_obj <- getHeatmapObject(ih)
#’   draw(hm_obj)
#’ }
#’ @export
setMethod(
  "InformativeHeatmap", "ANY",
  function(data,
           pch_val = 16,
           unit_val = 4,
           significant_color = "black",
           trending_color = "yellow",
           significant_pvalue = 0.05,
           trending_pvalue = 0.1,
           ...) 
  {
    ## 1. INPUT VALIDATION

    # Must be a numeric matrix
    if (!is.matrix(data)) {
      stop("`data` must be a matrix.")
    }
    if (!is.numeric(data)) {
      stop("`data` must contain only numeric values.")
    }
    # Dimensions
    n_rows <- nrow(data)
    n_cols <- ncol(data)
    if (n_rows < 1 || n_cols < 1) {
      stop("`data` must have at least one row and one column.")
    }

    # Validate overlay p‐value matrix (if provided via ...$significance_level)
    dots <- list(...)
    significance_level <- NULL
    if ("significance_level" %in% names(dots)) {
      significance_level <- dots$significance_level
      if (!is.matrix(significance_level) || !is.numeric(significance_level)) {
        stop("`significance_level` must be a numeric matrix.")
      }
      if (!identical(dim(significance_level), dim(data))) {
        stop(
          "`significance_level` dimensions (",
          paste(dim(significance_level), collapse = "×"),
          ") do not match `data` dimensions (",
          paste(dim(data), collapse = "×"),
          ")."
        )
      }
      # Remove it from dots so we don't pass it to Heatmap()
      dots$significance_level <- NULL
    }

    # Validate p‐value thresholds
    if (!is.numeric(significant_pvalue) || length(significant_pvalue) != 1 ||
        significant_pvalue <= 0 || significant_pvalue >= 1) {
      stop("`significant_pvalue` must be a single number in (0,1).")
    }
    if (!is.numeric(trending_pvalue) || length(trending_pvalue) != 1 ||
        trending_pvalue <= significant_pvalue || trending_pvalue >= 1) {
      stop("`trending_pvalue` must be a single number in (significant_pvalue, 1).")
    }

    # Validate pch_val and unit_val
    if (!is.numeric(pch_val) || length(pch_val) != 1) {
      stop("`pch_val` must be a single numeric or integer.")
    }
    if (!is.numeric(unit_val) || length(unit_val) != 1 || unit_val <= 0) {
      stop("`unit_val` must be a positive numeric (mm).")
    }

    # Check ComplexHeatmap availability
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop(
        "ComplexHeatmap is required. Install via BiocManager::install('ComplexHeatmap')."
      )
    }

    ## 2. BUILD CUSTOM LAYER FUNCTION (if p‐values provided)
    custom_layer_fun <- NULL
    if (!is.null(significance_level)) {
      # Pre‐cache locals for speed
      sig_p         <- significant_pvalue
      tr_p          <- trending_pvalue
      sig_col       <- significant_color
      tr_col        <- trending_color
      pmat_local    <- significance_level
      pch_local     <- pch_val
      unit_local    <- unit_val

      custom_layer_fun <- function(j, i, x, y, w, h, fill) {
        # j, i are integer vectors of column/row indices of each cell
        # restore_matrix returns a matrix of linear indices
        ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
        # ind_mat has same dims as the heatmap viewport; each element is linear index

        # Flatten once
        lin_idx <- as.vector(ind_mat)
        # Compute row/col indices
        # Given linear index k: row = ((k-1) %% n_rows) + 1, col = floor((k-1)/n_rows) + 1
        rows <- ((lin_idx - 1) %% n_rows) + 1
        cols <- floor((lin_idx - 1) / n_rows) + 1

        # Extract p‐values
        pv <- pmat_local[cbind(rows, cols)]

        # Determine colors: vectorized
        col_vec <- rep(NA_character_, length(pv))
        sig_idx <- which(pv < sig_p)
        tr_idx  <- which(pv >= sig_p & pv < tr_p)
        if (length(sig_idx) > 0) col_vec[sig_idx] <- sig_col
        if (length(tr_idx)  > 0) col_vec[tr_idx]  <- tr_col

        # Only draw points where col_vec is not NA
        keep <- which(!is.na(col_vec))
        if (length(keep) > 0) {
          grid::grid.points(
            x = x[lin_idx[keep]],
            y = y[lin_idx[keep]],
            pch = pch_local,
            gp = grid::gpar(col = col_vec[keep]),
            size = grid::unit(unit_local, "mm")
          )
        }
      }

      # Supply to Heatmap arguments
      dots$layer_fun <- custom_layer_fun
    }

    ## 3. CREATE COMPLEXHEATMAP OBJECT

    # Ensure data is passed as first argument
    hm_args <- c(list(matrix = data), dots)
    hm_obj  <- do.call(ComplexHeatmap::Heatmap, hm_args)

    ## 4. CONSTRUCT S4 OBJECT

    new(
      "InformativeHeatmap",
      heatmap = hm_obj,
      params  = list(
        pch_val             = pch_val,
        unit_val            = unit_val,
        significant_color   = significant_color,
        trending_color      = trending_color,
        significant_pvalue  = significant_pvalue,
        trending_pvalue     = trending_pvalue,
        extra_args          = dots
      )
    )
  }
)

#’ Define a Generic Method ‘updateLayerFun’
#’
#’ Update the overlay layer function of an existing InformativeHeatmap.
#’
#’ @param x An `InformativeHeatmap` object.
#’ @param layer_fun A function compatible with ComplexHeatmap’s layering system.
#’ @return Updated `InformativeHeatmap` object.
#’ @rdname updateLayerFun
#’ @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))

#’ @rdname updateLayerFun
#’ @export
setMethod(
  "updateLayerFun", "InformativeHeatmap",
  function(x, layer_fun) {
    if (!inherits(x, "InformativeHeatmap")) {
      stop("`x` must be an InformativeHeatmap object.")
    }
    if (!is.function(layer_fun)) {
      stop("`layer_fun` must be a function.")
    }
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop(
        "ComplexHeatmap is required. Install via BiocManager::install('ComplexHeatmap')."
      )
    }

    # Retrieve stored parameters
    old_args <- x@params$extra_args
    # Overwrite layer_fun
    old_args$layer_fun <- layer_fun

    # Recreate Heatmap with the same data
    mat_data <- x@heatmap@matrix
    hm_args   <- c(list(matrix = mat_data), old_args)
    new_hm    <- do.call(ComplexHeatmap::Heatmap, hm_args)

    x@heatmap <- new_hm
    x@params$extra_args <- old_args
    x
  }
)

#’ Create an InformativeHeatmap from a MultiAssayExperiment Object
#’
#’ Clusters samples using iClusterPlus, runs limma for each assay, combines logFCs,
#’ and overlays points colored by p‐value significance/trending.
#’
#’ @param mae A `MultiAssayExperiment` object containing ≥1 assay.
#’ @param significant_pvalue Numeric in (0,1) (default = 0.05).
#’ @param trending_pvalue Numeric in (significant_pvalue,1) (default = 0.1).
#’ @param significant_color Color for p < significant_pvalue (default = "black").
#’ @param trending_color Color for significant_pvalue ≤ p < trending_pvalue (default = "yellow").
#’ @param pch_val Integer or numeric pch for overlay (default = 16).
#’ @param unit_val Numeric point size in mm (default = 4).
#’ @param K Integer ≥2: number of clusters for iClusterPlus (default = 3).
#’ @param lambda Numeric in (0,∞): regularization for iClusterPlus (default = 0.2).
#’ @param coef Integer ≥1: which coefficient/contrast in limma to use (default = 2).
#’ @param ... Additional arguments passed to ComplexHeatmap::Heatmap.
#’
#’ @return An InformativeHeatmap object.
#’ @importFrom MultiAssayExperiment assays
#’ @importFrom stats model.matrix lmFit
#’ @importFrom limma eBayes topTable
#’ @export
InformativeHeatmapFromMAE <- function(
  mae,
  significant_pvalue = 0.05,
  trending_pvalue    = 0.1,
  significant_color  = "black",
  trending_color     = "yellow",
  pch_val            = 16,
  unit_val           = 4,
  K                  = 3,
  lambda             = 0.2,
  coef               = 2,
  ...
) {
  ## 1. INPUT VALIDATION

  # Check MultiAssayExperiment availability
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("MultiAssayExperiment package is required.")
  }
  # Check that mae is a MultiAssayExperiment
  if (!inherits(mae, "MultiAssayExperiment")) {
    stop("`mae` must be a MultiAssayExperiment.")
  }

  # Extract list of assay matrices (as numeric matrices)
  assay_list <- MultiAssayExperiment::assays(mae)
  if (length(assay_list) < 1) {
    stop("`mae` must contain at least one assay.")
  }
  # Convert each assay to numeric matrix and transpose (samples as rows)
  data_list <- lapply(assay_list, function(se) {
    mat <- as.matrix(se)
    if (!is.numeric(mat)) mode(mat) <- "numeric"
    t(mat)
  })
  # Name the list: dt1, dt2, …
  names(data_list) <- paste0("dt", seq_along(data_list))

  # Validate clustering parameters
  if (!requireNamespace("iClusterPlus", quietly = TRUE)) {
    stop("iClusterPlus is required. Install via BiocManager::install('iClusterPlus').")
  }
  if (!is.numeric(K) || length(K) != 1 || K < 2) {
    stop("`K` must be a single integer ≥ 2.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("`lambda` must be a single positive number.")
  }

  # Validate p-value thresholds
  if (!is.numeric(significant_pvalue) || length(significant_pvalue) != 1 ||
      significant_pvalue <= 0 || significant_pvalue >= 1) {
    stop("`significant_pvalue` must be a single number in (0,1).")
  }
  if (!is.numeric(trending_pvalue) || length(trending_pvalue) != 1 ||
      trending_pvalue <= significant_pvalue || trending_pvalue >= 1) {
    stop("`trending_pvalue` must be in (significant_pvalue, 1).")
  }
  if (!is.numeric(pch_val) || length(pch_val) != 1) {
    stop("`pch_val` must be a single numeric or integer.")
  }
  if (!is.numeric(unit_val) || length(unit_val) != 1 || unit_val <= 0) {
    stop("`unit_val` must be a positive numeric.")
  }
  if (!is.numeric(coef) || length(coef) != 1 || coef < 1) {
    stop("`coef` must be a single integer ≥ 1.")
  }

  ## 2. RUN iClusterPlus CLUSTERING

  # Prepare arguments
  icolist <- c(data_list, list(
    type   = rep("gaussian", length(data_list)),
    K      = as.integer(K),
    lambda = rep(lambda, length(data_list)),
    maxiter = as.integer(20)
  ))
  fit <- do.call(iClusterPlus::iClusterPlus, icolist)

  clusters <- factor(fit$clusters)
  if (length(clusters) != nrow(data_list[[1]])) {
    stop("Mismatch between number of samples in clusters and data.")
  }

  ## 3. DESIGN MATRIX FOR LIMMA

  design <- stats::model.matrix(~ 0 + clusters)
  colnames(design) <- levels(clusters)

  ## 4. RUN LIMMA FOR EACH ASSAY → STORE logFC AND p‐VALUES

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma is required. Install via BiocManager::install('limma').")
  }
  limma_results <- lapply(names(data_list), function(aname) {
    mat   <- data_list[[aname]]      # samples × features
    fitlm <- limma::lmFit(t(mat), design)  # transpose: features × samples
    fitlm <- limma::eBayes(fitlm)
    tt    <- limma::topTable(fitlm, coef = coef, number = Inf, sort.by = "none")
    if (!all(c("logFC", "P.Value") %in% colnames(tt))) {
      stop("topTable did not return logFC or P.Value.")
    }
    list(
      logFC    = tt$logFC,
      p_values = tt$P.Value
    )
  })
  names(limma_results) <- names(data_list)

  ## 5. COMBINE logFC & p‐VALUES ACROSS ASSAYS

  # Each element is a vector of length = number of rows in original assay
  logfc_list <- lapply(limma_results, function(res) res$logFC)
  pval_list  <- lapply(limma_results, function(res) res$p_values)

  # Combine by rbind → matrix assays × features
  combined_logFC    <- do.call(rbind, logfc_list)
  combined_pvalues  <- do.call(rbind, pval_list)

  if (!identical(dim(combined_logFC), dim(combined_pvalues))) {
    stop("Combined logFC and p-value matrices have mismatched dimensions.")
  }

  # Transpose → features × assays (rows × columns)
  heatmap_data     <- t(combined_logFC)
  pmat_for_overlay <- t(combined_pvalues)

  ## 6. BUILD CUSTOM LAYER FUNCTION FOR OVERLAY

  n_rows <- nrow(heatmap_data)
  n_cols <- ncol(heatmap_data)

  custom_layer_fun <- function(j, i, x, y, w, h, fill) {
    ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
    lin_idx <- as.vector(ind_mat)
    rows <- ((lin_idx - 1) %% n_rows) + 1
    cols <- floor((lin_idx - 1) / n_rows) + 1

    pv_vals <- pmat_for_overlay[cbind(rows, cols)]
    col_vec <- rep(NA_character_, length(pv_vals))
    sig_idx <- which(pv_vals < significant_pvalue)
    tr_idx  <- which(pv_vals >= significant_pvalue & pv_vals < trending_pvalue)
    if (length(sig_idx) > 0) col_vec[sig_idx] <- significant_color
    if (length(tr_idx) > 0) col_vec[tr_idx] <- trending_color

    keep <- which(!is.na(col_vec))
    if (length(keep) > 0) {
      grid::grid.points(
        x    = x[lin_idx[keep]],
        y    = y[lin_idx[keep]],
        pch  = pch_val,
        gp   = grid::gpar(col = col_vec[keep]),
        size = grid::unit(unit_val, "mm")
      )
    }
  }

  ## 7. CALL InformativeHeatmap CONSTRUCTOR

  ih <- InformativeHeatmap(
    heatmap_data,
    layer_fun           = custom_layer_fun,
    pch_val             = pch_val,
    unit_val            = unit_val,
    significant_color   = significant_color,
    trending_color      = trending_color,
    significant_pvalue  = significant_pvalue,
    trending_pvalue     = trending_pvalue,
    ...
  )
  return(ih)
}

#’ Generic to retrieve the underlying Heatmap object
#’
#’ @param x An InformativeHeatmap object.
#’ @return A ComplexHeatmap::Heatmap object.
#’ @rdname getHeatmapObject
#’ @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))

#’ Retrieve the Heatmap from InformativeHeatmap
#’
#’ @param x An InformativeHeatmap object.
#’ @return A ComplexHeatmap::Heatmap object.
#’ @export
setMethod(
  "getHeatmapObject", "InformativeHeatmap",
  function(x) {
    if (!inherits(x, "InformativeHeatmap")) {
      stop("`x` must be an InformativeHeatmap object.")
    }
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop(
        "ComplexHeatmap is required. Install via BiocManager::install('ComplexHeatmap')."
      )
    }
    return(x@heatmap)
  }
)






