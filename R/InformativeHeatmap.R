################################################################################
# InformativeHeatmap Class – Updated for Multimodal Integration
################################################################################

# Suppress notes about global variables
utils::globalVariables(c("assays", "setNames", "iClusterPlus"))

#’ -----------------------------------------------------------------------------
#’ InformativeHeatmap: A Class for Enhanced Heatmaps
#’
#’ Encapsulates a ComplexHeatmap::Heatmap object and its construction parameters.
#’ Now supports single‐modality (MAE, matrix, or ready‐to‐plot table) as well as
#’ multi‐modality (a named list of FC/p‐value matrices) in one integrative display.
#’
#’ @slot heatmap Heatmap object (or HeatmapList combined) from ComplexHeatmap
#’ @slot params  List of all parameters used (including sub‐structures for each modality)
#’ @exportClass InformativeHeatmap
#’ @import methods
#’ @importFrom ComplexHeatmap Heatmap
setClass(
  "InformativeHeatmap",
  slots = c(
    heatmap = "Heatmap",  # Could be a HeatmapList combined with ‘+’
    params  = "list"
  )
)

#’ -----------------------------------------------------------------------------
#’ Generic Constructor: InformativeHeatmap
#’ @param ... Arguments passed to specific methods
#’ @export
setGeneric("InformativeHeatmap", function(...) standardGeneric("InformativeHeatmap"))

#’ -----------------------------------------------------------------------------
#’ Constructor Method: Single‐Modality (Matrix + Metadata / MAE)
#’ @param data               Numeric matrix (features × samples). Or MultiAssayExperiment.
#’ @param meta               data.frame of sample metadata (rownames = colnames(data)).
#’                           Required if `data` is a matrix and `runClustering = FALSE`.
#’ @param runClustering      Logical: if TRUE, run iClusterPlus; else use `groupColumn`.
#’                           Default = FALSE.
#’ @param groupColumn        Character: column in `meta` for grouping (if runClustering = FALSE).
#’ @param continuous         One of "auto","continuous","count". Default = "auto".
#’ @param var_quantile       Numeric in [0,1] for variance filtering. NULL = no filtering.
#’ @param min_features       Integer: if var_quantile = NULL, keep top by variance. NULL = no filtering.
#’ @param pvalue_cutoff      Numeric (0–1). Default = 0.05.
#’ @param trending_cutoff    Numeric (0–1). Default = 0.1.
#’ @param fc_cutoff          Numeric ≥ 0. Minimum |log2FC| for DE. Default = 0.
#’ @param max_features       Integer. After DE filtering, keep at most this many features. NULL = no cap.
#’ @param significant_color  Color for p < pvalue_cutoff. Default = "red".
#’ @param trending_color     Color for p in [pvalue_cutoff, trending_cutoff). Default = "orange".
#’ @param pch_val            Integer. Plotting character. Default = 16.
#’ @param unit_val           Numeric. Point size in mm. Default = 2.
#’ @param K                  Integer. # clusters for iClusterPlus. Default = 3.
#’ @param lambda             Numeric or vector. Regularization for iClusterPlus. Default = 0.2.
#’ @param coef               Integer or character. Which coefficient in limma. Default = 2.
#’ @param BPPARAM            BiocParallelParam for parallel DE. Default = BiocParallel::bpparam().
#’ @param heatmap_data_scale "logFC" or "expression": basis for heatmap. Default = "logFC".
#’ @param cluster_rows       Logical. Default = TRUE.
#’ @param cluster_columns    Logical. Default = TRUE.
#’ @param show_row_names     Logical. Default = FALSE.
#’ @param show_column_names  Logical. Default = FALSE.
#’ @param col                Color mapping function. Default = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")).
#’ @param ...                Additional arguments passed to ComplexHeatmap::Heatmap().
#’ @return InformativeHeatmap S4 object
setMethod(
  "InformativeHeatmap",
  signature(data = "matrix"),
  function(
    data,
    meta                 = NULL,
    runClustering        = FALSE,
    groupColumn          = NULL,
    continuous           = c("auto","continuous","count"),
    var_quantile         = NULL,
    min_features         = NULL,
    pvalue_cutoff        = 0.05,
    trending_cutoff      = 0.1,
    fc_cutoff            = 0,
    max_features         = NULL,
    significant_color    = "red",
    trending_color       = "orange",
    pch_val              = 16,
    unit_val             = 2,
    K                    = 3,
    lambda               = 0.2,
    coef                 = 2,
    BPPARAM              = BiocParallel::bpparam(),
    heatmap_data_scale   = c("logFC","expression"),
    cluster_rows         = TRUE,
    cluster_columns      = TRUE,
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    col                  = circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red")),
    ...
  ) {
    ## 0) Argument Matching & Validation
    heatmap_data_scale <- match.arg(heatmap_data_scale)
    continuous         <- match.arg(continuous)
    if (!is.matrix(data) || !is.numeric(data)) {
      stop("`data` must be a numeric matrix (features × samples).")
    }
    if (any(is.na(data)) || any(is.infinite(data))) {
      stop("`data` contains NA or Inf. Remove or impute before calling.")
    }
    n_features <- nrow(data)
    n_samples  <- ncol(data)

    if (!runClustering) {
      if (is.null(meta) || !is.data.frame(meta)) {
        stop("When runClustering = FALSE, `meta` must be a data.frame of sample metadata.")
      }
      if (is.null(groupColumn) || !(groupColumn %in% colnames(meta))) {
        stop("`groupColumn` must be a valid column name in `meta` when runClustering = FALSE.")
      }
      if (!all(colnames(data) %in% rownames(meta))) {
        stop("Column names of `data` must exactly match row names of `meta`.")
      }
      meta <- meta[colnames(data), , drop = FALSE]
    }

    ## 1) Determine dataType
    if (continuous == "auto") {
      if (all(data == floor(data)) && max(data, na.rm = TRUE) > 30) {
        dataType <- "count"
      } else {
        dataType <- "continuous"
      }
    } else {
      dataType <- continuous
    }

    ## 2) Variance Filtering (optional)
    if (!is.null(var_quantile) && (var_quantile < 0 || var_quantile > 1)) {
      stop("`var_quantile` must be NULL or a numeric in [0,1].")
    }
    if (!is.null(var_quantile)) {
      vr <- matrixStats::rowVars(data, na.rm = TRUE)
      thr <- stats::quantile(vr, probs = var_quantile, na.rm = TRUE)
      keep_var <- vr >= thr
      if (!is.null(min_features) && sum(keep_var) < min_features) {
        ord_var <- order(vr, decreasing = TRUE)
        keep_var <- rep(FALSE, n_features)
        keep_var[ord_var[seq_len(min_features)]] <- TRUE
      }
      data_filt <- data[keep_var, , drop = FALSE]
      if (nrow(data_filt) == 0) {
        stop("No features remain after variance filtering—adjust `var_quantile` or `min_features`.")
      }
    } else if (!is.null(min_features)) {
      vr <- matrixStats::rowVars(data, na.rm = TRUE)
      ord_var <- order(vr, decreasing = TRUE)
      keep_var <- rep(FALSE, n_features)
      keep_var[ord_var[seq_len(min_features)]] <- TRUE
      data_filt <- data[keep_var, , drop = FALSE]
    } else {
      data_filt <- data
    }

    ## 3) Design Matrix: Clustering vs. Provided Group
    if (runClustering) {
      if (!requireNamespace("iClusterPlus", quietly = TRUE)) {
        stop("iClusterPlus is required for clustering; install it via BiocManager::install('iClusterPlus').")
      }
      dt    <- t(data_filt)  # samples × features
      icl   <- tryCatch({
        iClusterPlus::iClusterPlus(
          dt1     = dt,
          type    = "gaussian",
          K       = K,
          lambda  = rep(lambda, 1),
          maxiter = 20
        )
      }, error = function(e) {
        stop("iClusterPlus error: ", e$message)
      })
      clusters <- factor(icl$clusters)
      if (length(levels(clusters)) < 2) {
        stop("iClusterPlus produced fewer than 2 clusters; adjust `K` or filtering.")
      }
      design <- stats::model.matrix(~ 0 + clusters)
      colnames(design) <- levels(clusters)
    } else {
      grp <- meta[[groupColumn]]
      if (!is.factor(grp)) grp <- factor(grp)
      if (length(levels(grp)) < 2) {
        stop("`groupColumn` must have at least two levels.")
      }
      design <- stats::model.matrix(~ 0 + grp)
      colnames(design) <- levels(grp)
    }

    ## 4) Differential Expression (limma or voom+limma)
    if (dataType == "count") {
      v   <- limma::voom(data_filt, design, plot = FALSE)
      fit <- limma::lmFit(v, design)
    } else {
      fit <- limma::lmFit(data_filt, design)
    }
    fit <- limma::eBayes(fit)

    ## 5) Identify Coefficient
    if (is.character(coef)) {
      ci <- which(colnames(design) == coef)
      if (length(ci) != 1) {
        stop("`coef` (character) did not match any design column.")
      }
      coef_idx <- ci
    } else if (is.numeric(coef)) {
      coef_idx <- as.integer(coef)
      if (coef_idx < 1 || coef_idx > ncol(design)) {
        stop("`coef` (numeric) out of range for design.")
      }
    } else {
      stop("`coef` must be integer or character.")
    }

    ## 6) topTable & Basic Filtering
    tt <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
    keep_valid <- !is.na(tt$logFC) & !is.na(tt$P.Value) & is.finite(tt$logFC) & is.finite(tt$P.Value)
    tt <- tt[keep_valid, , drop = FALSE]
    if (nrow(tt) == 0) {
      stop("No valid DE results after filtering NA/Inf in topTable.")
    }

    ## 7) Apply fc_cutoff
    if (fc_cutoff > 0) {
      keep_fc <- abs(tt$logFC) >= fc_cutoff
      tt       <- tt[keep_fc, , drop = FALSE]
      if (nrow(tt) == 0) {
        stop("No features pass the `fc_cutoff` threshold.")
      }
    }

    ## 8) Apply max_features by |logFC|
    if (!is.null(max_features) && max_features < nrow(tt)) {
      ord_fc <- order(abs(tt$logFC), decreasing = TRUE)
      tt     <- tt[ord_fc[seq_len(max_features)], , drop = FALSE]
    }

    ## 9) Prepare Heatmap Data & P-Value Vector
    if (heatmap_data_scale == "logFC") {
      heatmap_mat <- matrix(tt$logFC, nrow = 1,
                            dimnames = list("logFC", rownames(tt)))
      pvals_vec   <- tt$P.Value
    } else {  # "expression"
      sel_feats   <- rownames(tt)
      heatmap_mat <- data_filt[sel_feats, , drop = FALSE]
      pvals_vec   <- tt$P.Value
    }

    ## 10) Vectorized layer_fun
    layer_fun <- function(j, i, x, y, w, h, fill) {
      if (heatmap_data_scale == "logFC") {
        # Single‐row heatmap: i == 1 for all; j indexes columns/features
        pv_vec   <- pvals_vec[j]
        cols_all <- ifelse(pv_vec < pvalue_cutoff, significant_color,
                    ifelse(pv_vec < trending_cutoff, trending_color, NA))
        keep_idx <- which(!is.na(cols_all))
        if (length(keep_idx) > 0) {
          grid::grid.points(
            x    = x[keep_idx],
            y    = y[keep_idx],
            pch  = pch_val,
            gp   = grid::gpar(col = cols_all[keep_idx]),
            size = grid::unit(unit_val, "mm")
          )
        }
      } else {
        # Multi-row heatmap: i indexes feature‐rows, j indexes samples
        pv_all   <- pvals_vec[i]
        cols_all <- ifelse(pv_all < pvalue_cutoff, significant_color,
                    ifelse(pv_all < trending_cutoff, trending_color, NA))
        keep_idx <- which(!is.na(cols_all))
        if (length(keep_idx) > 0) {
          grid::grid.points(
            x    = x[keep_idx],
            y    = y[keep_idx],
            pch  = pch_val,
            gp   = grid::gpar(col = cols_all[keep_idx]),
            size = grid::unit(unit_val, "mm")
          )
        }
      }
    }

    ## 11) Build the Heatmap
    ht <- ComplexHeatmap::Heatmap(
      heatmap_mat,
      name              = "Value",
      col               = col,
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
      show_column_names = show_column_names,
      layer_fun         = layer_fun,
      ...
    )

    ## 12) Package parameters
    params_list <- list(
      data               = data,
      meta               = if (!runClustering) meta else NULL,
      runClustering      = runClustering,
      groupColumn        = if (!runClustering) groupColumn else NULL,
      continuous         = dataType,
      var_quantile       = var_quantile,
      min_features       = min_features,
      pvalue_cutoff      = pvalue_cutoff,
      trending_cutoff    = trending_cutoff,
      fc_cutoff          = fc_cutoff,
      max_features       = max_features,
      significant_color  = significant_color,
      trending_color     = trending_color,
      pch_val            = pch_val,
      unit_val           = unit_val,
      K                  = if (runClustering) K else NA,
      lambda             = if (runClustering) lambda else NA,
      coef               = coef,
      BPPARAM            = BPPARAM,
      heatmap_data_scale = heatmap_data_scale,
      cluster_rows       = cluster_rows,
      cluster_columns    = cluster_columns,
      show_row_names     = show_row_names,
      show_column_names  = show_column_names,
      col                = col
    )

    methods::new(
      "InformativeHeatmap",
      heatmap = ht,
      params  = params_list
    )
  }
)

#’ -----------------------------------------------------------------------------
#’ Constructor Method: MAE Workflow
#’ Delegates to the matrix method after extracting assay + metadata.
#’ @param data                MultiAssayExperiment
#’ @param assayName          Character or integer: which assay to use
#’ @param ...                Passed to InformativeHeatmap(matrix, ...)
#’ @return InformativeHeatmap S4 object
setMethod(
  "InformativeHeatmap",
  signature(data = "MultiAssayExperiment"),
  function(
    data,
    assayName            = NULL,
    runClustering        = FALSE,
    groupColumn          = NULL,
    continuous           = c("auto","continuous","count"),
    var_quantile         = NULL,
    min_features         = NULL,
    pvalue_cutoff        = 0.05,
    trending_cutoff      = 0.1,
    fc_cutoff            = 0,
    max_features         = NULL,
    significant_color    = "red",
    trending_color       = "orange",
    pch_val              = 16,
    unit_val             = 2,
    K                    = 3,
    lambda               = 0.2,
    coef                 = 2,
    BPPARAM              = BiocParallel::bpparam(),
    heatmap_data_scale   = c("logFC","expression"),
    cluster_rows         = TRUE,
    cluster_columns      = TRUE,
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    col                  = circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red")),
    ...
  ) {
    continuous <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)

    if (!inherits(data, "MultiAssayExperiment")) {
      stop("`data` must be a MultiAssayExperiment when using this method.")
    }
    assays_list <- MultiAssayExperiment::experiments(data)
    if (is.null(assayName)) {
      if (length(assays_list) != 1) {
        stop("MAE has multiple assays; please specify `assayName`.")
      }
      assayName <- names(assays_list)[1]
    }
    if (!(assayName %in% names(assays_list))) {
      stop("Assay '", assayName, "' not found in MAE.")
    }
    se   <- assays_list[[assayName]]
    expr <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)

    # Delegate to matrix‐based constructor
    callGeneric(
      data               = expr,
      meta               = meta,
      runClustering      = runClustering,
      groupColumn        = groupColumn,
      continuous         = continuous,
      var_quantile       = var_quantile,
      min_features       = min_features,
      pvalue_cutoff      = pvalue_cutoff,
      trending_cutoff    = trending_cutoff,
      fc_cutoff          = fc_cutoff,
      max_features       = max_features,
      significant_color  = significant_color,
      trending_color     = trending_color,
      pch_val            = pch_val,
      unit_val           = unit_val,
      K                  = K,
      lambda             = lambda,
      coef               = coef,
      BPPARAM            = BPPARAM,
      heatmap_data_scale = heatmap_data_scale,
      cluster_rows       = cluster_rows,
      cluster_columns    = cluster_columns,
      show_row_names     = show_row_names,
      show_column_names  = show_column_names,
      col                = col,
      ...
    )
  }
)

#’ -----------------------------------------------------------------------------
#’ Constructor: Ready‐to‐Plot Table Workflow
#’ @param fc_matrix          Numeric matrix (features × contrasts) of log2FC.
#’ @param pval_matrix        Numeric matrix (same dims) of p‐values.
#’ @param pvalue_cutoff      Numeric: p < this → “significant.” Default = 0.05.
#’ @param trending_cutoff    Numeric: p ∈ [pvalue_cutoff, trending_cutoff) → “trending.” Default = 0.1.
#’ @param fc_cutoff          Numeric: keep only rows where max|fc| ≥ fc_cutoff. Default = 0.
#’ @param max_features       Integer: after filtering, keep top by max|fc|. NULL = no cap.
#’ @param significant_color  Color for p < pvalue_cutoff. Default = "red".
#’ @param trending_color     Color for p ∈ [pvalue_cutoff, trending_cutoff). Default = "orange".
#’ @param pch_val            Integer: plotting character. Default = 16.
#’ @param unit_val           Numeric: point size (mm). Default = 2.
#’ @param cluster_rows       Logical. Default = TRUE.
#’ @param cluster_columns    Logical. Default = TRUE.
#’ @param show_row_names     Logical. Default = TRUE.
#’ @param show_column_names  Logical. Default = TRUE.
#’ @param col                Color mapping function. Default = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")).
#’ @param ...                Additional arguments for ComplexHeatmap::Heatmap().
#’ @return InformativeHeatmap S4 object
setGeneric("InformativeHeatmap_table", function(fc_matrix, pval_matrix, ...) standardGeneric("InformativeHeatmap_table"))

setMethod(
  "InformativeHeatmap_table",
  signature(fc_matrix = "matrix", pval_matrix = "matrix"),
  function(
    fc_matrix,
    pval_matrix,
    pvalue_cutoff      = 0.05,
    trending_cutoff    = 0.1,
    fc_cutoff          = 0,
    max_features       = NULL,
    significant_color  = "red",
    trending_color     = "orange",
    pch_val            = 16,
    unit_val           = 2,
    cluster_rows       = TRUE,
    cluster_columns    = TRUE,
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    col                = circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red")),
    ...
  ) {
    ## 0) Validations
    if (!all(dim(fc_matrix) == dim(pval_matrix))) {
      stop("`fc_matrix` and `pval_matrix` must have identical dimensions.")
    }
    if (!is.numeric(fc_matrix) || !is.numeric(pval_matrix)) {
      stop("`fc_matrix` and `pval_matrix` must be numeric matrices.")
    }
    if (any(is.na(fc_matrix)) || any(is.infinite(fc_matrix))) {
      stop("`fc_matrix` contains NA or Inf; remove or impute first.")
    }
    if (any(is.na(pval_matrix)) || any(is.infinite(pval_matrix))) {
      stop("`pval_matrix` contains NA or Inf; remove or impute first.")
    }

    n_feat <- nrow(fc_matrix)

    ## 1) Filter by fc_cutoff
    if (fc_cutoff > 0) {
      max_fc_per_row <- matrixStats::rowMaxs(abs(fc_matrix))
      keep_fc         <- max_fc_per_row >= fc_cutoff
      if (sum(keep_fc) == 0) {
        stop("No features pass the `fc_cutoff` threshold.")
      }
      fc_matrix   <- fc_matrix[keep_fc, , drop = FALSE]
      pval_matrix <- pval_matrix[keep_fc, , drop = FALSE]
    }

    ## 2) Keep top max_features by max |FC|
    if (!is.null(max_features) && max_features < nrow(fc_matrix)) {
      max_fc_per_row <- matrixStats::rowMaxs(abs(fc_matrix))
      ord_idx        <- order(max_fc_per_row, decreasing = TRUE)
      keep_idx       <- ord_idx[seq_len(max_features)]
      fc_matrix      <- fc_matrix[keep_idx, , drop = FALSE]
      pval_matrix    <- pval_matrix[keep_idx, , drop = FALSE]
    }

    ## 3) Vectorized layer_fun
    layer_fun_tab <- function(j, i, x, y, w, h, fill) {
      # i = row indices, j = column indices
      pv_vec   <- pval_matrix[cbind(i, j)]
      cols_all <- rep(NA_character_, length(pv_vec))
      sig_idx  <- which(pv_vec < pvalue_cutoff)
      trend_idx<- which(pv_vec >= pvalue_cutoff & pv_vec < trending_cutoff)
      cols_all[sig_idx]   <- significant_color
      cols_all[trend_idx] <- trending_color
      keep_idx <- which(!is.na(cols_all))
      if (length(keep_idx) > 0) {
        grid::grid.points(
          x    = x[keep_idx],
          y    = y[keep_idx],
          pch  = pch_val,
          gp   = grid::gpar(col = cols_all[keep_idx]),
          size = grid::unit(unit_val, "mm")
        )
      }
    }

    ## 4) Build Heatmap
    ht_tab <- ComplexHeatmap::Heatmap(
      fc_matrix,
      name              = "logFC",
      col               = col,
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
      show_column_names = show_column_names,
      layer_fun         = layer_fun_tab,
      ...
    )

    params_list <- list(
      fc_matrix         = fc_matrix,
      pval_matrix       = pval_matrix,
      pvalue_cutoff     = pvalue_cutoff,
      trending_cutoff   = trending_cutoff,
      fc_cutoff         = fc_cutoff,
      max_features      = max_features,
      significant_color = significant_color,
      trending_color    = trending_color,
      pch_val           = pch_val,
      unit_val          = unit_val,
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
      show_column_names = show_column_names,
      col               = col
    )

    methods::new(
      "InformativeHeatmap",
      heatmap = ht_tab,
      params  = params_list
    )
  }
)

#’ -----------------------------------------------------------------------------
#’ Constructor: Multimodal Integration (List of FC/P‐Value Matrices)
#’
#’ Users supply a named list of FC matrices and a corresponding named list of
#’ p‐value matrices. Each modality yields its own heatmap + overlay; all are
#’ combined side-by-side into a single ComplexHeatmap “+” object.
#’
#’ @param fc_list            Named list of numeric matrices (features × samples) of log2FC.
#’ @param pval_list          Named list of numeric matrices (same dims as fc_list) of p‐values.
#’ @param pvalue_cutoff      Numeric: p < this → “significant.” Default = 0.05.
#’ @param trending_cutoff    Numeric: p ∈ [pvalue_cutoff, trending_cutoff) → “trending.” Default = 0.1.
#’ @param fc_cutoff          Numeric ≥0: keep only rows where max|fc| ≥ fc_cutoff. Default = 0.
#’ @param max_features       Integer: after filtering each modality, keep top by max|fc|. NULL = no cap.
#’ @param significant_color  Color for p < pvalue_cutoff. Default = "red".
#’ @param trending_color     Color for p ∈ [pvalue_cutoff, trending_cutoff). Default = "orange".
#’ @param pch_val            Integer: plotting character for overlay. Default = 16.
#’ @param unit_val           Numeric: point size in mm. Default = 2.
#’ @param cluster_rows       Logical. Cluster rows. Default = TRUE.
#’ @param cluster_columns    Logical. Cluster columns. Default = TRUE.
#’ @param show_row_names     Logical. Default = TRUE.
#’ @param show_column_names  Logical. Default = TRUE.
#’ @param col_list           Named list of color mapping functions, same names as fc_list.
#’ @param ...                Additional arguments for ComplexHeatmap::Heatmap() for each modality.
#’ @return InformativeHeatmap S4 object combining all modalities.
#’ @examples
#’ \dontrun{
#’ fm <- list(
#’   RNA  = fc_mat_rna,
#’   Methyl = fc_mat_meth
#’ )
#’ pm <- list(
#’   RNA  = pval_mat_rna,
#’   Methyl = pval_mat_meth
#’ )
#’ ih_multi <- InformativeHeatmap(
#’   fc_list          = fm,
#’   pval_list        = pm,
#’   pvalue_cutoff    = 0.01,
#’   trending_cutoff  = 0.05,
#’   fc_cutoff        = 1,
#’   max_features     = 100,
#’   significant_color = "darkred",
#’   trending_color    = "gold",
#’   pch_val           = 16,
#’   unit_val          = 2,
#’   cluster_rows      = TRUE,
#’   cluster_columns   = TRUE,
#’   show_row_names    = TRUE,
#’   show_column_names = TRUE,
#’   col_list          = list(
#’     RNA     = circlize::colorRamp2(c(-3,0,3), c("navy","white","firebrick")),
#’     Methyl  = circlize::colorRamp2(c(-2,0,2), c("green","white","purple"))
#’   )
#’ )
#’ ComplexHeatmap::draw(getHeatmapObject(ih_multi))
#’ }
setMethod(
  "InformativeHeatmap",
  signature(data = "list"),
  function(
    fc_list,
    pval_list,
    pvalue_cutoff      = 0.05,
    trending_cutoff    = 0.1,
    fc_cutoff          = 0,
    max_features       = NULL,
    significant_color  = "red",
    trending_color     = "orange",
    pch_val            = 16,
    unit_val           = 2,
    cluster_rows       = TRUE,
    cluster_columns    = TRUE,
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    col_list           = NULL,
    ...
  ) {
    ## 0) Validate Input Lists
    if (!is.list(fc_list) || !is.list(pval_list)) {
      stop("`fc_list` and `pval_list` must both be named lists of matrices.")
    }
    if (!identical(sort(names(fc_list)), sort(names(pval_list)))) {
      stop("Names of `fc_list` and `pval_list` must match exactly.")
    }
    modality_names <- names(fc_list)
    if (is.null(modality_names) || any(modality_names == "")) {
      stop("`fc_list` must be a named list (each element has a nonempty name).")
    }

    # For each modality, produce one Heatmap
    hm_list <- list()
    subparams_list <- list()

    for (mod in modality_names) {
      fc_mat  <- fc_list[[mod]]
      pv_mat  <- pval_list[[mod]]

      # Validate each pair
      if (!is.matrix(fc_mat) || !is.numeric(fc_mat)) {
        stop(sprintf("fc_list[['%s']] must be a numeric matrix.", mod))
      }
      if (!is.matrix(pv_mat) || !is.numeric(pv_mat)) {
        stop(sprintf("pval_list[['%s']] must be a numeric matrix.", mod))
      }
      if (!all(dim(fc_mat) == dim(pv_mat))) {
        stop(sprintf("Dimensions of fc_list[['%s']] and pval_list[['%s']] must match.", mod, mod))
      }
      if (any(is.na(fc_mat)) || any(is.infinite(fc_mat))) {
        stop(sprintf("fc_list[['%s']] contains NA or Inf; remove or impute first.", mod))
      }
      if (any(is.na(pv_mat)) || any(is.infinite(pv_mat))) {
        stop(sprintf("pval_list[['%s']] contains NA or Inf; remove or impute first.", mod))
      }

      # 1) fc_cutoff filtering
      nr <- nrow(fc_mat)
      if (fc_cutoff > 0) {
        row_max_fc <- matrixStats::rowMaxs(abs(fc_mat))
        keep_row   <- row_max_fc >= fc_cutoff
        if (sum(keep_row) == 0) {
          stop(sprintf("No features in modality '%s' pass fc_cutoff=%.2f", mod, fc_cutoff))
        }
        fc_mat  <- fc_mat[keep_row, , drop = FALSE]
        pv_mat  <- pv_mat[keep_row, , drop = FALSE]
      }

      # 2) max_features filtering
      if (!is.null(max_features) && max_features < nrow(fc_mat)) {
        row_max_fc <- matrixStats::rowMaxs(abs(fc_mat))
        ord_idx    <- order(row_max_fc, decreasing = TRUE)
        keep_idx   <- ord_idx[seq_len(max_features)]
        fc_mat     <- fc_mat[keep_idx, , drop = FALSE]
        pv_mat     <- pv_mat[keep_idx, , drop = FALSE]
      }

      # 3) Determine color function for this modality
      if (!is.null(col_list) && mod %in% names(col_list)) {
        col_fun_mod <- col_list[[mod]]
      } else {
        col_fun_mod <- circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red"))
      }

      # 4) Build layer_fun for this modality
      layer_fun_mod <- function(j, i, x, y, w, h, fill) {
        pv_vec   <- pv_mat[cbind(i, j)]
        cols_all <- rep(NA_character_, length(pv_vec))
        sig_idx  <- which(pv_vec < pvalue_cutoff)
        tr_idx   <- which(pv_vec >= pvalue_cutoff & pv_vec < trending_cutoff)
        cols_all[sig_idx] <- significant_color
        cols_all[tr_idx]  <- trending_color
        keep_idx <- which(!is.na(cols_all))
        if (length(keep_idx) > 0) {
          grid::grid.points(
            x    = x[keep_idx],
            y    = y[keep_idx],
            pch  = pch_val,
            gp   = grid::gpar(col = cols_all[keep_idx]),
            size = grid::unit(unit_val, "mm")
          )
        }
      }

      # 5) Create the Heatmap for this modality
      ht_mod <- ComplexHeatmap::Heatmap(
        fc_mat,
        name              = paste0(mod, "_logFC"),
        col               = col_fun_mod,
        cluster_rows      = cluster_rows,
        cluster_columns   = cluster_columns,
        show_row_names    = show_row_names,
        show_column_names = show_column_names,
        layer_fun         = layer_fun_mod,
        ...
      )

      hm_list[[mod]] <- ht_mod
      subparams_list[[mod]] <- list(
        fc_matrix         = fc_mat,
        pval_matrix       = pv_mat,
        pvalue_cutoff     = pvalue_cutoff,
        trending_cutoff   = trending_cutoff,
        fc_cutoff         = fc_cutoff,
        max_features      = max_features,
        significant_color = significant_color,
        trending_color    = trending_color,
        pch_val           = pch_val,
        unit_val          = unit_val,
        cluster_rows      = cluster_rows,
        cluster_columns   = cluster_columns,
        show_row_names    = show_row_names,
        show_column_names = show_column_names,
        col_fun           = col_fun_mod
      )
    }

    # 6) Combine all modality‐specific Heatmaps using '+'
    combined_ht <- Reduce(function(a, b) a + b, hm_list)

    ## 7) Package parameters
    params_list <- list(
      multimodal        = TRUE,
      modalities        = modality_names,
      subparams         = subparams_list
    )

    methods::new(
      "InformativeHeatmap",
      heatmap = combined_ht,
      params  = params_list
    )
  }
)

#’ -----------------------------------------------------------------------------
#’ Retrieve the underlying Heatmap (or HeatmapList) from an InformativeHeatmap
#’ @param x InformativeHeatmap object
#’ @return Heatmap or combined HeatmapList
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))

#’ @rdname getHeatmapObject
#’ @export
setMethod(
  "getHeatmapObject",
  signature(x = "InformativeHeatmap"),
  function(x) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop("ComplexHeatmap is required to retrieve the Heatmap object.")
    }
    return(x@heatmap)
  }
)

#’ -----------------------------------------------------------------------------
#’ Update the layer_fun of an existing InformativeHeatmap
#’ @param x         InformativeHeatmap object
#’ @param layer_fun New layer_fun (compatible with ComplexHeatmap)
#’ @return Updated InformativeHeatmap
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))

#’ @rdname updateLayerFun
#’ @export
setMethod(
  "updateLayerFun",
  signature(x = "InformativeHeatmap"),
  function(x, layer_fun) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop("ComplexHeatmap is required to update layer function.")
    }
    params <- x@params

    # If multimodal, we cannot generically update across all modalities;
    # user should reconstruct with new layer_fun for each modality manually.
    if (isTRUE(params$multimodal)) {
      stop("To update layer_fun for a multimodal InformativeHeatmap, please reconstruct modality by modality.")
    }

    # Single‐modality case: retrieve stored parameters, replace layer_fun, rebuild
    hm_mat <- x@heatmap@matrix
    hm_args <- c(list(hm_mat), params)
    hm_args$layer_fun <- layer_fun
    new_ht <- do.call(ComplexHeatmap::Heatmap, hm_args)

    x@heatmap <- new_ht
    x@params  <- modifyList(params, list(layer_fun = layer_fun))
    return(x)
  }
)

