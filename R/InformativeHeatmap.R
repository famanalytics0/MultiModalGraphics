################################################################################
# InformativeHeatmap Class – Updated for Multimodal Integration
################################################################################

# Suppress notes about global variables
utils::globalVariables(c("assays_list", "iClusterPlus"))

#’ @title InformativeHeatmap: A Class for Enhanced Heatmaps
#’ @description
#’ Encapsulates a ComplexHeatmap::Heatmap (or combined HeatmapList) and the
#’ parameters used to build it, supporting single-modality (matrix/MAE/table)
#’ and multimodal (named lists of FC/p-value matrices) workflows.
#’
#’ @slot heatmap A ComplexHeatmap::Heatmap or combined HeatmapList
#’ @slot params  A list of all parameters passed (including sub-lists for each modality)
#’ @exportClass InformativeHeatmap
#’ @import methods
#’ @importFrom ComplexHeatmap Heatmap
setClass(
  "InformativeHeatmap",
  slots = c(
    heatmap = "Heatmap",
    params  = "list"
  )
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @param data Either a numeric matrix (features × samples), a MultiAssayExperiment,
#’             or a named list of FC matrices.
#’ @param ...  Further args passed to the matching constructor method.
#’ @export
setGeneric(
  "InformativeHeatmap",
  function(data, ...) standardGeneric("InformativeHeatmap")
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @describeIn InformativeHeatmap Single-modality constructor from matrix + metadata
#’ @param meta               data.frame of sample metadata (rownames = colnames(data))
#’ @param runClustering      Logical; run iClusterPlus if TRUE
#’ @param groupColumn        Column in `meta` to define groups if !runClustering
#’ @param continuous         One of "auto","continuous","count"
#’ @param var_quantile       Quantile filter for variance
#’ @param min_features       If no quantile filtering, keep top N by variance
#’ @param pvalue_cutoff      p-value threshold for “significance”
#’ @param trending_cutoff    p-value threshold for “trending”
#’ @param fc_cutoff          Minimum |log2FC| for DE
#’ @param max_features       After DE filter, cap features
#’ @param significant_color  Color for significant p
#’ @param trending_color     Color for trending p
#’ @param pch_val            Point character
#’ @param unit_val           Point size (mm)
#’ @param K                  # clusters for iClusterPlus
#’ @param lambda             Regularization for iClusterPlus
#’ @param coef               Which coefficient in limma
#’ @param BPPARAM            BiocParallelParam
#’ @param heatmap_data_scale One of "logFC","expression"
#’ @param cluster_rows       Logical
#’ @param cluster_columns    Logical
#’ @param show_row_names     Logical
#’ @param show_column_names  Logical
#’ @param col                Color mapping function (from circlize)
#’ @param ...                Other args to ComplexHeatmap::Heatmap
#’ @export
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
    col                  = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    ...
  ) {
    ## --- Argument checks & matching ---
    continuous         <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)
    if (!is.matrix(data) || !is.numeric(data))
      stop("`data` must be a numeric matrix.")
    if (anyNA(data) || any(is.infinite(data)))
      stop("`data` contains NA or Inf.")
    n_features <- nrow(data)
    ## Metadata for grouping
    if (!runClustering) {
      if (!is.data.frame(meta))
        stop("`meta` must be provided as a data.frame when runClustering = FALSE.")
      if (!groupColumn %in% colnames(meta))
        stop("`groupColumn` not found in `meta`.")
      if (!all(colnames(data) %in% rownames(meta)))
        stop("Column names of `data` must match row names of `meta`.")
      meta <- meta[colnames(data), , drop = FALSE]
    }
    ## --- Variance filtering ---
    if (!is.null(var_quantile)) {
      if (var_quantile < 0 || var_quantile > 1)
        stop("`var_quantile` must be in [0,1].")
      vr <- matrixStats::rowVars(data, na.rm = TRUE)
      thr <- stats::quantile(vr, var_quantile, na.rm = TRUE)
      keep <- vr >= thr
      if (!is.null(min_features) && sum(keep) < min_features) {
        ord <- order(vr, decreasing = TRUE)
        keep[] <- FALSE; keep[ord[seq_len(min_features)]] <- TRUE
      }
      data <- data[keep, , drop = FALSE]
      if (!nrow(data)) stop("No features remain after variance filtering.")
    } else if (!is.null(min_features)) {
      vr  <- matrixStats::rowVars(data, na.rm = TRUE)
      ord <- order(vr, decreasing = TRUE)
      data <- data[ord[seq_len(min_features)], , drop = FALSE]
    }
    ## --- Design matrix ---
    if (runClustering) {
      if (!requireNamespace("iClusterPlus", quietly=TRUE))
        stop("Install 'iClusterPlus' for clustering.")
      ic <- iClusterPlus::iClusterPlus(
        dt1     = t(data),
        type    = "gaussian",
        K       = K,
        lambda  = rep(lambda, length = length(lambda)),
        maxiter = 20
      )
      clusters <- factor(ic$clusters)
      design <- stats::model.matrix(~ 0 + clusters)
      colnames(design) <- levels(clusters)
    } else {
      grp <- meta[[groupColumn]]
      if (!is.factor(grp)) grp <- factor(grp)
      if (nlevels(grp) < 2)
        stop("`groupColumn` must have ≥2 levels.")
      design <- stats::model.matrix(~ 0 + grp)
      colnames(design) <- levels(grp)
    }
    ## --- Differential Expression ---
    if (continuous=="count") {
      vfit <- limma::voom(data, design, plot=FALSE)
      fit  <- limma::lmFit(vfit, design)
    } else {
      fit  <- limma::lmFit(data, design)
    }
    fit <- limma::eBayes(fit)
    ## Determine `coef_idx`
    if (is.character(coef)) {
      coef_idx <- match(coef, colnames(design))
      if (is.na(coef_idx)) stop("`coef` not found in design.")
    } else {
      coef_idx <- as.integer(coef)
      if (coef_idx<1 || coef_idx>ncol(design))
        stop("`coef` index out of range.")
    }
    tt <- limma::topTable(fit, coef=coef_idx, number=Inf, sort.by="none")
    keep_valid <- complete.cases(tt[,c("logFC","P.Value")])
    tt <- tt[keep_valid, , drop=FALSE]
    if (!nrow(tt)) stop("No valid DE results.")
    if (fc_cutoff>0) {
      keep_fc <- abs(tt$logFC) >= fc_cutoff
      tt <- tt[keep_fc, , drop=FALSE]
      if (!nrow(tt)) stop("No features pass `fc_cutoff`.")
    }
    if (!is.null(max_features) && max_features < nrow(tt)) {
      ord <- order(abs(tt$logFC), decreasing=TRUE)
      tt  <- tt[ord[seq_len(max_features)], , drop=FALSE]
    }
    ## --- Prepare heatmap matrix + overlay points ---
    if (heatmap_data_scale=="logFC") {
      mat    <- matrix(tt$logFC, nrow=1,
                       dimnames=list("logFC", rownames(tt)))
      pvals  <- tt$P.Value
    } else {
      feats  <- rownames(tt)
      mat    <- data[feats, , drop=FALSE]
      pvals  <- tt$P.Value
    }
    layer_fun <- function(j, i, x, y, w, h, fill) {
      if (heatmap_data_scale=="logFC") {
        pv_vec <- pvals[j]
        cols   <- ifelse(pv_vec < pvalue_cutoff, significant_color,
                         ifelse(pv_vec < trending_cutoff, trending_color, NA))
        idx    <- which(!is.na(cols))
        if (length(idx))
          grid::grid.points(x=x[idx], y=y[idx], pch=pch_val,
                            gp=grid::gpar(col=cols[idx]),
                            size=grid::unit(unit_val,"mm"))
      } else {
        pv_vec <- pvals[i]
        cols   <- ifelse(pv_vec < pvalue_cutoff, significant_color,
                         ifelse(pv_vec < trending_cutoff, trending_color, NA))
        idx    <- which(!is.na(cols))
        if (length(idx))
          grid::grid.points(x=x[idx], y=y[idx], pch=pch_val,
                            gp=grid::gpar(col=cols[idx]),
                            size=grid::unit(unit_val,"mm"))
      }
    }
    ## --- Build heatmap ---
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
    ## --- Collect params & return ---
    params <- as.list(environment())[c(
      "runClustering","groupColumn","continuous","var_quantile","min_features",
      "pvalue_cutoff","trending_cutoff","fc_cutoff","max_features",
      "significant_color","trending_color","pch_val","unit_val",
      "K","lambda","coef","BPPARAM","heatmap_data_scale",
      "cluster_rows","cluster_columns","show_row_names",
      "show_column_names","col"
    )]
    methods::new("InformativeHeatmap", heatmap=ht, params=params)
  }
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @describeIn InformativeHeatmap MAE constructor that extracts assay & metadata
#’ @param data       A `MultiAssayExperiment`
#’ @param assayName  Name or index of assay to use; auto-selects if only one
#’ @export
setMethod(
  "InformativeHeatmap",
  signature(data = "MultiAssayExperiment"),
  function(
    data,
    assayName            = NULL,
    ...
  ) {
    if (!inherits(data, "MultiAssayExperiment"))
      stop("Input must be a MultiAssayExperiment.")
    assays_list <- MultiAssayExperiment::experiments(data)
    if (is.null(assayName)) {
      if (length(assays_list)!=1)
        stop("Multiple assays found; please specify `assayName`.")
      assayName <- names(assays_list)[1]
    }
    if (!assayName %in% names(assays_list))
      stop("Assay '", assayName, "' not in MAE.")
    se   <- assays_list[[assayName]]
    expr <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors=FALSE)
    callGeneric(data=expr, meta=meta, ...)
  }
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @describeIn InformativeHeatmap Multimodal constructor from named lists
#’ @param data      Named list of log2FC matrices
#’ @param pval_list Named list of matching p-value matrices
#’ @export
setMethod(
  "InformativeHeatmap",
  signature(data = "list"),
  function(
    data,
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
    ## Validate lists
    if (!is.list(data) || !is.list(pval_list))
      stop("`data` and `pval_list` must be named lists.")
    if (!identical(sort(names(data)), sort(names(pval_list))))
      stop("Names of FC and p-value lists must match.")
    mods <- names(data)
    hm_list <- vector("list", length(mods)); names(hm_list) <- mods
    params_sub <- vector("list", length(mods)); names(params_sub) <- mods

    for (mod in mods) {
      fc  <- data[[mod]]; pv <- pval_list[[mod]]
      ## same-dimension and numeric checks
      if (!is.matrix(fc) || !is.numeric(fc)) stop("FC for ",mod," not numeric matrix.")
      if (!is.matrix(pv) || !is.numeric(pv)) stop("p-vals for ",mod," not numeric matrix.")
      if (!all(dim(fc)==dim(pv)))
        stop("Dimensions of FC/p-val for ",mod," differ.")
      ## filter by fc_cutoff, max_features (vectorized)
      nr <- nrow(fc)
      if (fc_cutoff>0) {
        mfc <- matrixStats::rowMaxs(abs(fc))
        keep <- mfc>=fc_cutoff
        fc   <- fc[keep,,drop=FALSE]; pv <- pv[keep,,drop=FALSE]
        if (!nrow(fc)) stop("No features pass fc_cutoff in ",mod)
      }
      if (!is.null(max_features) && max_features < nrow(fc)) {
        mfc  <- matrixStats::rowMaxs(abs(fc))
        ord  <- order(mfc, decreasing=TRUE)[1:max_features]
        fc   <- fc[ord,,drop=FALSE]; pv <- pv[ord,,drop=FALSE]
      }
      ## color function
      col_fun <- if (!is.null(col_list) && mod%in%names(col_list))
                   col_list[[mod]]
                 else circlize::colorRamp2(c(-2,0,2), c("blue","white","red"))
      ## build layer_fun
      lf_mod <- function(j,i,x,y,...) {
        pvij   <- pv[cbind(i,j)]
        cols   <- rep(NA_character_, length(pvij))
        cols[pvij<pvalue_cutoff] <- significant_color
        idx2   <- which(pvij>=pvalue_cutoff & pvij<trending_cutoff)
        cols[idx2] <- trending_color
        idx    <- which(!is.na(cols))
        if (length(idx))
          grid::grid.points(x=x[idx], y=y[idx], pch=pch_val,
                            gp=grid::gpar(col=cols[idx]),
                            size=grid::unit(unit_val,"mm"))
      }
      hm_list[[mod]] <- ComplexHeatmap::Heatmap(
        fc,
        name              = paste0(mod,"_logFC"),
        col               = col_fun,
        cluster_rows      = cluster_rows,
        cluster_columns   = cluster_columns,
        show_row_names    = show_row_names,
        show_column_names = show_column_names,
        layer_fun         = lf_mod,
        ...
      )
      params_sub[[mod]] <- list(fc=fc,pv=pv)
    }
    combined <- Reduce(`+`, hm_list)
    params <- list(
      multimodal        = TRUE,
      subparams         = params_sub,
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
      show_column_names = show_column_names
    )
    methods::new("InformativeHeatmap", heatmap=combined, params=params)
  }
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @export
setGeneric(
  "InformativeHeatmap_table",
  function(fc_matrix, pval_matrix, ...) standardGeneric("InformativeHeatmap_table")
)
#’ @rdname InformativeHeatmap
#' @export
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
    col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    ...
  ) {
    if (!all(dim(fc_matrix)==dim(pval_matrix)))
      stop("Dimensions of fc_matrix and pval_matrix must match.")
    if (!is.numeric(fc_matrix)||!is.numeric(pval_matrix))
      stop("Both matrices must be numeric.")
    ## filtering by fc_cutoff & max_features...
    nr <- nrow(fc_matrix)
    if (fc_cutoff>0) {
      mfc  <- matrixStats::rowMaxs(abs(fc_matrix))
      keep <- mfc>=fc_cutoff
      fc_matrix <- fc_matrix[keep,,drop=FALSE]
      pval_matrix <- pval_matrix[keep,,drop=FALSE]
      if (!nrow(fc_matrix)) stop("No features pass fc_cutoff.")
    }
    if (!is.null(max_features) && max_features < nrow(fc_matrix)) {
      mfc  <- matrixStats::rowMaxs(abs(fc_matrix))
      ord  <- order(mfc, decreasing=TRUE)[1:max_features]
      fc_matrix <- fc_matrix[ord,,drop=FALSE]
      pval_matrix <- pval_matrix[ord,,drop=FALSE]
    }
    ## build layer_fun
    layer_fun_tab <- function(j,i,x,y,...) {
      pv   <- pval_matrix[cbind(i,j)]
      cols <- rep(NA_character_, length(pv))
      cols[pv<pvalue_cutoff] <- significant_color
      idx2 <- which(pv>=pvalue_cutoff & pv<trending_cutoff)
      cols[idx2] <- trending_color
      idx <- which(!is.na(cols))
      if (length(idx))
        grid::grid.points(x=x[idx], y=y[idx], pch=pch_val,
                          gp=grid::gpar(col=cols[idx]),
                          size=grid::unit(unit_val,"mm"))
    }
    ht <- ComplexHeatmap::Heatmap(
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
    params <- as.list(environment())[c(
      "pvalue_cutoff","trending_cutoff","fc_cutoff","max_features",
      "significant_color","trending_color","pch_val","unit_val",
      "cluster_rows","cluster_columns","show_row_names","show_column_names","col"
    )]
    methods::new("InformativeHeatmap", heatmap=ht, params=params)
  }
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))
#’ @rdname InformativeHeatmap
setMethod(
  "getHeatmapObject",
  signature(x = "InformativeHeatmap"),
  function(x) x@heatmap
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))
#’ @rdname InformativeHeatmap
setMethod(
  "updateLayerFun",
  signature(x = "InformativeHeatmap"),
  function(x, layer_fun) {
    params <- x@params
    if (isTRUE(params$multimodal))
      stop("Rebuild multimodal heatmap to change layer_fun.")
    args <- c(list(x@heatmap@matrix), params)
    args$layer_fun <- layer_fun
    new_ht <- do.call(ComplexHeatmap::Heatmap, args)
    x@heatmap <- new_ht
    x@params  <- modifyList(params, list(layer_fun=layer_fun))
    x
  }
)

#’ -----------------------------------------------------------------------------
#’ @rdname InformativeHeatmap
#’ @exportMethod show
setMethod(
  "show",
  signature(object = "InformativeHeatmap"),
  function(object) {
    ComplexHeatmap::draw(object@heatmap)
    invisible(object)
  }
)

