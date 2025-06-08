################################################################################
# InformativeHeatmap.R
################################################################################

# Suppress notes about global variables
utils::globalVariables(c("assays_list", "iClusterPlus"))

#' -----------------------------------------------------------------------------
#' @title InformativeHeatmap: A Class for Enhanced Heatmaps
#' @description
#' Encapsulates a ComplexHeatmap::Heatmap (or combined HeatmapList) and the
#' parameters used to build it. Supports:
#'  - single-modality from a matrix + metadata
#'  - MAE workflows
#'  - ready-to-plot FC/p-value tables
#'  - multimodal integration (named lists of FC/p-value matrices)
#'
#' @slot heatmap A ComplexHeatmap::Heatmap or HeatmapList
#' @slot params  A list of all parameters used (including sub-lists for each modality)
#' @exportClass InformativeHeatmap
#' @importMethodsFrom methods setClass setGeneric setMethod new
#' @importFrom ComplexHeatmap Heatmap
setClass(
  "InformativeHeatmap",
  slots = c(
    heatmap = "Heatmap",
    params  = "list"
  )
)

#' -----------------------------------------------------------------------------
#' @rdname InformativeHeatmap
#' @param data Either
#'   * a numeric matrix (features × samples),
#'   * a MultiAssayExperiment,
#'   * a named list of FC matrices (for multimodal),
#'   * or a matrix + p-value matrix (for table workflow via InformativeHeatmap_table())
#' @param ... Passed to the matching method
#' @export
setGeneric(
  "InformativeHeatmap",
  function(data, ...) standardGeneric("InformativeHeatmap")
)

#' -----------------------------------------------------------------------------
#' @describeIn InformativeHeatmap Single-modality constructor from matrix + metadata
#' @param meta               data.frame of sample metadata (rownames = colnames(data))
#' @param runClustering      Logical; run iClusterPlus if TRUE
#' @param groupColumn        Column in `meta` to define groups if !runClustering
#' @param continuous         One of "auto","continuous","count"
#' @param var_quantile       Quantile filter for variance (NULL = no filter)
#' @param min_features       If no quantile filter, keep top N by variance (NULL = no cap)
#' @param pvalue_cutoff      p-value threshold for significance
#' @param trending_cutoff    p-value threshold for trending
#' @param fc_cutoff          Minimum |log2FC| for DE
#' @param max_features       After DE, cap by top |FC| (NULL = no cap)
#' @param significant_color  Color for p < pvalue_cutoff
#' @param trending_color     Color for p in [pvalue_cutoff, trending_cutoff)
#' @param pch_val            Point character (integer)
#' @param unit_val           Point size (mm)
#' @param K                  # clusters for iClusterPlus
#' @param lambda             Regularization for iClusterPlus
#' @param coef               Which coefficient in limma (integer or name)
#' @param BPPARAM            BiocParallelParam
#' @param heatmap_data_scale One of "logFC","expression"
#' @param cluster_rows       Logical
#' @param cluster_columns    Logical
#' @param show_row_names     Logical
#' @param show_column_names  Logical
#' @param col                Color mapping function (e.g. circlize::colorRamp2)
#' @param ...                Other args for ComplexHeatmap::Heatmap
#' @export
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
    ## --- argument matching & validation ---
    continuous         <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)
    if (!is.matrix(data) || !is.numeric(data))
      stop("`data` must be a numeric matrix.")
    if (anyNA(data) || any(is.infinite(data)))
      stop("`data` contains NA or Inf.")
    nfeat <- nrow(data)

    ## metadata checks
    if (!runClustering) {
      if (!is.data.frame(meta))
        stop("`meta` must be a data.frame when runClustering = FALSE.")
      if (!groupColumn %in% colnames(meta))
        stop("`groupColumn` not found in `meta`.")
      if (!all(colnames(data) %in% rownames(meta)))
        stop("Column names of `data` must match row names of `meta`.")
      meta <- meta[colnames(data), , drop = FALSE]
    }

    ## variance filtering
    if (!is.null(var_quantile)) {
      if (var_quantile < 0 || var_quantile > 1)
        stop("`var_quantile` must be in [0,1].")
      vr  <- matrixStats::rowVars(data, na.rm = TRUE)
      thr <- stats::quantile(vr, probs = var_quantile, na.rm = TRUE)
      keep <- vr >= thr
      if (!is.null(min_features) && sum(keep) < min_features) {
        ord <- order(vr, decreasing = TRUE)
        keep[] <- FALSE
        keep[ord[seq_len(min_features)]] <- TRUE
      }
      data <- data[keep, , drop = FALSE]
      if (nrow(data)==0) stop("No features remain after variance filtering.")
    } else if (!is.null(min_features)) {
      vr  <- matrixStats::rowVars(data, na.rm = TRUE)
      ord <- order(vr, decreasing = TRUE)
      data <- data[ord[seq_len(min_features)], , drop = FALSE]
    }

    ## design matrix
    if (runClustering) {
      if (!requireNamespace("iClusterPlus", quietly=TRUE))
        stop("Install 'iClusterPlus' for clustering.")
      ic <- iClusterPlus::iClusterPlus(
        dt1     = t(data),
        type    = "gaussian",
        K       = K,
        lambda  = rep(lambda, length=lambda),
        maxiter = 20
      )
      clusters <- factor(ic$clusters)
      design <- stats::model.matrix(~0+clusters)
      colnames(design) <- levels(clusters)
    } else {
      grp <- meta[[groupColumn]]
      if (!is.factor(grp)) grp <- factor(grp)
      if (nlevels(grp) < 2) stop("`groupColumn` must have ≥2 levels.")
      design <- stats::model.matrix(~0+grp)
      colnames(design) <- levels(grp)
    }

    ## differential expression
    fit <- if (continuous=="count") {
      vfit <- limma::voom(data, design, plot=FALSE)
      limma::lmFit(vfit, design)
    } else {
      limma::lmFit(data, design)
    }
    fit <- limma::eBayes(fit)

    ## choose coefficient
    coef_idx <- if (is.character(coef)) {
      which(colnames(design)==coef) |> 
        (\(x) if (length(x)!=1) stop("`coef` not found."))()
    } else {
      as.integer(coef)
    }
    tt <- limma::topTable(fit, coef=coef_idx, number=Inf, sort.by="none")
    keep_valid <- complete.cases(tt[,c("logFC","P.Value")])
    tt <- tt[keep_valid, , drop=FALSE]
    if (nrow(tt)==0) stop("No valid DE results.")
    if (fc_cutoff>0) {
      tt <- tt[abs(tt$logFC) >= fc_cutoff, , drop=FALSE]
      if (nrow(tt)==0) stop("No features pass `fc_cutoff`.")
    }
    if (!is.null(max_features) && nrow(tt)>max_features) {
      ord <- order(abs(tt$logFC), decreasing=TRUE)
      tt  <- tt[ord[seq_len(max_features)], , drop=FALSE]
    }

    ## prepare heatmap matrix + overlay p-values
    if (heatmap_data_scale=="logFC") {
      mat   <- matrix(tt$logFC, nrow=1,
                      dimnames=list("logFC", rownames(tt)))
      pvals <- tt$P.Value
    } else {
      feats <- rownames(tt)
      mat   <- data[feats, , drop=FALSE]
      pvals <- tt$P.Value
    }

    layer_fun <- function(j, i, x, y, w, h, fill) {
      pv_vec <- if (heatmap_data_scale=="logFC") pvals[j] else pvals[i]
      cols   <- ifelse(pv_vec < pvalue_cutoff, significant_color,
                       ifelse(pv_vec < trending_cutoff, trending_color, NA))
      idx    <- which(!is.na(cols))
      if (length(idx))
        grid::grid.points(
          x    = x[idx],
          y    = y[idx],
          pch  = pch_val,
          gp   = grid::gpar(col = cols[idx]),
          size = grid::unit(unit_val, "mm")
        )
    }

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name              = "Value",
      col               = col,
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
      show_column_names = show_column_names,
      layer_fun         = layer_fun,
      ...
    )

    params <- as.list(environment())[setdiff(names(formals()), c("data","..."))]
    methods::new("InformativeHeatmap", heatmap=ht, params=params)
  }
)

#' -----------------------------------------------------------------------------
#' @describeIn InformativeHeatmap MAE constructor — delegates to matrix method
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "MultiAssayExperiment"),
  function(data, assayName=NULL, ...) {
    assays_list <- MultiAssayExperiment::experiments(data)
    if (is.null(assayName)) {
      if (length(assays_list)!=1)
        stop("Multiple assays: specify `assayName`.")
      assayName <- names(assays_list)[1]
    }
    if (!(assayName %in% names(assays_list)))
      stop("Assay '", assayName, "' not found.")
    se   <- assays_list[[assayName]]
    expr <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors=FALSE)
    callGeneric(expr, meta=meta, ...)
  }
)

#' -----------------------------------------------------------------------------
#' @describeIn InformativeHeatmap Multimodal integration from named lists
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "list"),
  function(data, pval_list, ...) {
    if (!is.list(pval_list))
      stop("`pval_list` must be a list matching `data`.")
    if (!identical(sort(names(data)), sort(names(pval_list))))
      stop("Names of FC and p-value lists must match.")
    hm_list <- lapply(names(data), function(mod) {
      fc <- data[[mod]]
      pv <- pval_list[[mod]]
      if (!is.matrix(fc) || !is.matrix(pv))
        stop("Each element must be a numeric matrix.")
      # delegate to table constructor internally:
      InformativeHeatmap_table(
        fc_matrix        = fc,
        pval_matrix      = pv,
        name             = paste0(mod, "_logFC"),
        ...
      )
    })
    # extract the Heatmap objects and combine
    combined <- Reduce(function(a,b) a + b, lapply(hm_list, slot, "heatmap"))
    params   <- list(multimodal=TRUE, modalities=names(data))
    methods::new("InformativeHeatmap", heatmap=combined, params=params)
  }
)

#' -----------------------------------------------------------------------------
#' @title InformativeHeatmap_table: Table-based constructor
#' @describeIn InformativeHeatmap Build from an FC matrix + p-value matrix
#' @param name             Name for the heatmap legend (default = "logFC")
#' @export
setGeneric(
  "InformativeHeatmap_table",
  function(fc_matrix, pval_matrix, ...) standardGeneric("InformativeHeatmap_table")
)

#' @rdname InformativeHeatmap_table
#' @export
setMethod(
  "InformativeHeatmap_table",
  signature(fc_matrix="matrix", pval_matrix="matrix"),
  function(
    fc_matrix,
    pval_matrix,
    name               = "logFC",
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
    col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    ...
  ) {
    # input checks
    if (!all(dim(fc_matrix)==dim(pval_matrix)))
      stop("fc_matrix and pval_matrix must have identical dimensions.")
    if (!is.numeric(fc_matrix)||!is.numeric(pval_matrix))
      stop("Both matrices must be numeric.")
    # fc_cutoff filter
    if (fc_cutoff>0) {
      mfc <- matrixStats::rowMaxs(abs(fc_matrix))
      keep <- mfc >= fc_cutoff
      fc_matrix   <- fc_matrix[keep,,drop=FALSE]
      pval_matrix <- pval_matrix[keep,,drop=FALSE]
      if (nrow(fc_matrix)==0) stop("No features pass fc_cutoff.")
    }
    # max_features cap
    if (!is.null(max_features) && max_features < nrow(fc_matrix)) {
      mfc <- matrixStats::rowMaxs(abs(fc_matrix))
      ord <- order(mfc, decreasing=TRUE)[seq_len(max_features)]
      fc_matrix   <- fc_matrix[ord,,drop=FALSE]
      pval_matrix <- pval_matrix[ord,,drop=FALSE]
    }
    # build layer_fun
    layer_fun <- function(j, i, x, y, w, h, fill) {
      pv_vec <- pval_matrix[cbind(i,j)]
      cols   <- rep(NA_character_, length(pv_vec))
      cols[pv_vec < pvalue_cutoff] <- significant_color
      idx2   <- which(pv_vec >= pvalue_cutoff & pv_vec < trending_cutoff)
      cols[idx2] <- trending_color
      idx   <- which(!is.na(cols))
      if (length(idx))
        grid::grid.points(
          x    = x[idx], y = y[idx],
          pch  = pch_val,
          gp   = grid::gpar(col=cols[idx]),
          size = grid::unit(unit_val,"mm")
        )
    }
    ht <- ComplexHeatmap::Heatmap(
      fc_matrix,
      name              = name,
      col               = col,
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
      show_column_names = show_column_names,
      layer_fun         = layer_fun,
      ...
    )
    params <- as.list(environment())[setdiff(names(formals()), c("fc_matrix","pval_matrix","..."))]
    methods::new("InformativeHeatmap", heatmap=ht, params=params)
  }
)

#' -----------------------------------------------------------------------------
#' @rdname InformativeHeatmap
#' @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))

#' @rdname InformativeHeatmap
#' @export
setMethod("getHeatmapObject", "InformativeHeatmap", function(x) x@heatmap)

#' -----------------------------------------------------------------------------
#' @rdname InformativeHeatmap
#' @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))

#' @rdname InformativeHeatmap
#' @export
setMethod("updateLayerFun", "InformativeHeatmap", function(x, layer_fun) {
  if (isTRUE(x@params$multimodal))
    stop("Rebuild multimodal heatmap to change layer_fun.")
  args <- c(list(x@heatmap@matrix), x@params)
  args$layer_fun <- layer_fun
  new_ht <- do.call(ComplexHeatmap::Heatmap, args)
  x@heatmap <- new_ht
  x@params  <- modifyList(x@params, list(layer_fun=layer_fun))
  x
})

#' -----------------------------------------------------------------------------
#' @rdname InformativeHeatmap
#' @exportMethod show
setMethod("show", "InformativeHeatmap", function(object) {
  ComplexHeatmap::draw(object@heatmap)
  invisible(object)
})

################################################################################
#' @title InformativeHeatmapFromMAE: Convenient Wrapper for MAE Inputs
#' @description
#' Extracts an assay and its sample metadata from a `MultiAssayExperiment`, then
#' delegates to the matrix+metadata constructor of `InformativeHeatmap()`.
#'
#' @param mae       A `MultiAssayExperiment`
#' @param assayName Character or integer; which assay to use (auto-select if one)
#' @param ...       Passed on to `InformativeHeatmap(matrix, meta, ...)`
#' @return An object of class `InformativeHeatmap`
#' @export
InformativeHeatmapFromMAE <- function(mae, assayName=NULL, ...) {
  if (!inherits(mae, "MultiAssayExperiment"))
    stop("`mae` must be a MultiAssayExperiment.")
  assays_list <- MultiAssayExperiment::experiments(mae)
  if (is.null(assayName)) {
    if (length(assays_list)!=1) stop("Multiple assays: specify `assayName`.")
    assayName <- names(assays_list)[1]
  }
  if (!(assayName %in% names(assays_list)))
    stop("Assay '", assayName, "' not found.")
  se   <- assays_list[[assayName]]
  expr <- SummarizedExperiment::assay(se)
  meta <- as.data.frame(SummarizedExperiment::colData(se),
                        stringsAsFactors=FALSE)
  InformativeHeatmap(expr, meta=meta, ...)
}
