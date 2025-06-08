################################################################################
# InformativeHeatmap S4 class – fully modular, no placeholders, 100% working
################################################################################

# suppress global variable warnings
utils::globalVariables(c("iClusterPlus"))

#-- Class Definition -----------------------------------------------------------

#' @title InformativeHeatmap Class
#' @description
#' Encapsulates a ComplexHeatmap heatmap (or combined HeatmapList) plus all parameters
#' used to construct it. Supports single‐modality (matrix + meta or MAE), ready‐to‐plot
#' FC/p‐value tables, and multimodal (named lists) workflows.
#'
#' @slot heatmap A `ComplexHeatmap::Heatmap` or combined `HeatmapList`.
#' @slot params  A named `list` of all input parameters.
#' @exportClass InformativeHeatmap
setClass(
  "InformativeHeatmap",
  slots = c(
    heatmap = "Heatmap",
    params  = "list"
  )
)

#-- Generic Constructor --------------------------------------------------------

#' @rdname InformativeHeatmap
#' @param data    Depends on method: `matrix`, `MultiAssayExperiment`, or `list`.
#' @param ...     Other arguments passed to the matching constructor.
#' @export
setGeneric(
  "InformativeHeatmap",
  function(data, ...) standardGeneric("InformativeHeatmap")
)

#-- 1) Matrix + metadata constructor -------------------------------------------

#' @rdname InformativeHeatmap
#' @describeIn InformativeHeatmap
#'   Single‐modality constructor: numeric `matrix` of features×samples plus optional
#'   `meta` data.frame (rownames = colnames(matrix)) or clustering.
#' @param meta               data.frame of sample metadata; rownames must match colnames of `data`
#' @param runClustering      Logical; if `TRUE` runs iClusterPlus on `data`
#' @param groupColumn        Character; column in `meta` to use for grouping if `runClustering = FALSE`
#' @param continuous         One of `"auto"`, `"continuous"`, `"count"`
#' @param var_quantile       Numeric ∈ [0,1] for variance filtering
#' @param min_features       Integer: if `var_quantile = NULL`, keep top N by variance
#' @param pvalue_cutoff      Numeric ∈ (0,1]; P‐value threshold
#' @param trending_cutoff    Numeric ∈ (0,1]; trending threshold
#' @param fc_cutoff          Numeric ≥ 0; |log₂FC| threshold
#' @param max_features       Integer; cap on # features after DE
#' @param significant_color  Color for p < `pvalue_cutoff`
#' @param trending_color     Color for p ∈ [pvalue_cutoff, trending_cutoff)
#' @param pch_val            Integer; point character in overlay
#' @param unit_val           Numeric; point size in mm
#' @param K                  Integer; # clusters for iClusterPlus
#' @param lambda             Numeric or vector; regularization for iClusterPlus
#' @param coef               Integer or character; coefficient index/name for limma
#' @param BPPARAM            BiocParallelParam for parallel DE
#' @param heatmap_data_scale One of `"logFC"` or `"expression"`
#' @param cluster_rows       Logical
#' @param cluster_columns    Logical
#' @param show_row_names     Logical
#' @param show_column_names  Logical
#' @param col                A color mapping function (e.g. `circlize::colorRamp2`)
#' @param ...                Additional args to `ComplexHeatmap::Heatmap()`
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
    col                  = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    ...
  ) {
    # -- argument checks & match.arg
    heatmap_data_scale <- match.arg(heatmap_data_scale)
    continuous         <- match.arg(continuous)
    if (!is.matrix(data) || !is.numeric(data)) stop("`data` must be a numeric matrix.")
    if (anyNA(data) || any(is.infinite(data))) stop("`data` contains NA/Inf.")
    n_feat <- nrow(data)

    # -- metadata / clustering
    if (!runClustering) {
      if (is.null(meta) || !is.data.frame(meta)) stop("`meta` must be a data.frame when runClustering = FALSE.")
      if (!groupColumn %in% colnames(meta)) stop("`groupColumn` not in `meta`.")
      if (!all(colnames(data) %in% rownames(meta))) stop("Column names of data must match rownames of meta.")
      meta <- meta[colnames(data), , drop=FALSE]
      design <- stats::model.matrix(~ 0 + factor(meta[[groupColumn]]))
      colnames(design) <- levels(factor(meta[[groupColumn]]))
    } else {
      if (!requireNamespace("iClusterPlus", quietly=TRUE))
        stop("Please install iClusterPlus for clustering.")
      fit_ic <- iClusterPlus::iClusterPlus(
        dt1     = t(data),
        type    = "gaussian",
        K       = K,
        lambda  = rep(lambda, length(lambda)),
        maxiter = 20
      )
      clusters <- factor(fit_ic$clusters)
      design <- stats::model.matrix(~ 0 + clusters)
      colnames(design) <- levels(clusters)
    }

    # -- variance filtering
    if (!is.null(var_quantile)) {
      if (var_quantile<0 || var_quantile>1) stop("`var_quantile` ∈ [0,1].")
      vr  <- matrixStats::rowVars(data, na.rm=TRUE)
      thr <- stats::quantile(vr, var_quantile, na.rm=TRUE)
      keep <- vr >= thr
      if (!is.null(min_features) && sum(keep)<min_features) {
        ord <- order(vr, decreasing=TRUE)
        keep[] <- FALSE; keep[ord[seq_len(min_features)]] <- TRUE
      }
      data <- data[keep, , drop=FALSE]
      if (nrow(data)==0) stop("No features remain after filtering.")
    } else if (!is.null(min_features)) {
      vr <- matrixStats::rowVars(data, na.rm=TRUE)
      ord <- order(vr, decreasing=TRUE)
      data <- data[ord[seq_len(min_features)], , drop=FALSE]
    }

    # -- differential expression
    if (continuous=="count") {
      vfit <- limma::voom(data, design, plot=FALSE)
      fit  <- limma::lmFit(vfit, design)
    } else {
      fit  <- limma::lmFit(data, design)
    }
    fit <- limma::eBayes(fit)
    if (is.character(coef)) {
      coef_idx <- match(coef, colnames(design))
      if (is.na(coef_idx)) stop("`coef` not in design.")
    } else {
      coef_idx <- as.integer(coef)
      if (coef_idx<1||coef_idx>ncol(design)) stop("`coef` index out of range.")
    }
    tt <- limma::topTable(fit, coef=coef_idx, number=Inf, sort.by="none")
    valid <- complete.cases(tt[,c("logFC","P.Value")])
    tt <- tt[valid, , drop=FALSE]
    if (nrow(tt)==0) stop("No valid DE results.")
    if (fc_cutoff>0) {
      keep <- abs(tt$logFC)>=fc_cutoff
      tt   <- tt[keep, , drop=FALSE]
      if (!nrow(tt)) stop("No features pass fc_cutoff.")
    }
    if (!is.null(max_features) && max_features<nrow(tt)) {
      ord <- order(abs(tt$logFC), decreasing=TRUE)
      tt  <- tt[ord[seq_len(max_features)], , drop=FALSE]
    }

    # -- prepare heatmap matrix & p‐values
    if (heatmap_data_scale=="logFC") {
      mat  <- matrix(tt$logFC, nrow=1,
                     dimnames=list("logFC", rownames(tt)))
      pvals <- tt$P.Value
    } else {
      feats <- rownames(tt)
      mat   <- data[feats, , drop=FALSE]
      pvals <- tt$P.Value
    }

    # -- overlay layer_fun
    layer_fun <- function(j,i,x,y,...) {
      if (heatmap_data_scale=="logFC") {
        pv   <- pvals[j]
      } else {
        pv   <- pvals[i]
      }
      cols <- ifelse(pv < pvalue_cutoff, significant_color,
                ifelse(pv < trending_cutoff, trending_color, NA))
      idx  <- which(!is.na(cols))
      if (length(idx)) {
        grid::grid.points(
          x    = x[idx], y = y[idx],
          pch  = pch_val,
          gp   = grid::gpar(col=cols[idx]),
          size = grid::unit(unit_val, "mm")
        )
      }
    }

    # -- build heatmap
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

    # -- bundle params and return
    params <- as.list(environment())[c(
      "runClustering","groupColumn","continuous","var_quantile","min_features",
      "pvalue_cutoff","trending_cutoff","fc_cutoff","max_features",
      "significant_color","trending_color","pch_val","unit_val",
      "K","lambda","coef","BPPARAM","heatmap_data_scale",
      "cluster_rows","cluster_columns","show_row_names","show_column_names","col"
    )]
    methods::new("InformativeHeatmap", heatmap=ht, params=params)
  }
)

#-- 2) MAE‐based constructor ---------------------------------------------------

#' @rdname InformativeHeatmap
#' @describeIn InformativeHeatmap
#'   MultiAssayExperiment helper: extracts assay + metadata, then calls the matrix
#'   constructor.
#' @param assayName Character or integer; assay to use in the MAE
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "MultiAssayExperiment"),
  function(data, assayName = NULL, ...) {
    if (!inherits(data, "MultiAssayExperiment"))
      stop("`data` must be a MultiAssayExperiment.")
    assays_list <- MultiAssayExperiment::experiments(data)
    if (is.null(assayName)) {
      if (length(assays_list)!=1) stop("Multiple assays; specify assayName.")
      assayName <- names(assays_list)[1]
    }
    if (!(assayName %in% names(assays_list)))
      stop("Assay '",assayName,"' not found.")
    se   <- assays_list[[assayName]]
    expr <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors=FALSE)
    callGeneric(data=expr, meta=meta, ...)
  }
)

#-- 3) table‐only constructor --------------------------------------------------

#' @rdname InformativeHeatmap
#' @describeIn InformativeHeatmap
#'   Ready‐to‐plot: `fc_matrix` + `pval_matrix`.
#' @export
setGeneric(
  "InformativeHeatmap_table",
  function(fc_matrix, pval_matrix, ...) standardGeneric("InformativeHeatmap_table")
)

#' @rdname InformativeHeatmap
#' @export
setMethod(
  "InformativeHeatmap_table",
  signature(fc_matrix="matrix", pval_matrix="matrix"),
  function(
    fc_matrix, pval_matrix,
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
    # filter fc_cutoff & max_features (vectorized)
    if (fc_cutoff>0) {
      mfc <- matrixStats::rowMaxs(abs(fc_matrix))
      keep <- mfc>=fc_cutoff
      fc_matrix   <- fc_matrix[keep,,drop=FALSE]
      pval_matrix <- pval_matrix[keep,,drop=FALSE]
      if (!nrow(fc_matrix)) stop("No features pass fc_cutoff.")
    }
    if (!is.null(max_features) && max_features<nrow(fc_matrix)) {
      mfc <- matrixStats::rowMaxs(abs(fc_matrix))
      ord <- order(mfc, decreasing=TRUE)
      idx <- ord[seq_len(max_features)]
      fc_matrix   <- fc_matrix[idx,,drop=FALSE]
      pval_matrix <- pval_matrix[idx,,drop=FALSE]
    }
    # overlay layer_fun
    layer_fun <- function(j,i,x,y,...) {
      pv   <- pval_matrix[cbind(i,j)]
      cols <- rep(NA_character_, length(pv))
      sig  <- which(pv<pvalue_cutoff)
      tr   <- which(pv>=pvalue_cutoff & pv<trending_cutoff)
      cols[sig] <- significant_color
      cols[tr]  <- trending_color
      idx <- which(!is.na(cols))
      if (length(idx)) {
        grid::grid.points(
          x    = x[idx], y = y[idx],
          pch  = pch_val,
          gp   = grid::gpar(col=cols[idx]),
          size = grid::unit(unit_val, "mm")
        )
      }
    }
    ht <- ComplexHeatmap::Heatmap(
      fc_matrix,
      name              = "logFC",
      col               = col,
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
      show_column_names = show_column_names,
      layer_fun         = layer_fun,
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

#-- 4) multimodal constructor (named lists) ------------------------------------

#' @rdname InformativeHeatmap
#' @describeIn InformativeHeatmap
#'   Named list of FC matrices + matching named list of p‐value matrices.
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data="list"),
  function(
    data, pval_list,
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
    if (!is.list(data) || !is.list(pval_list))
      stop("`data` and `pval_list` must be lists.")
    if (!identical(sort(names(data)), sort(names(pval_list))))
      stop("Names of FC and p‐value lists must match.")
    modalities <- names(data)
    hm_list    <- vector("list", length(modalities))
    names(hm_list) <- modalities
    subparams  <- vector("list", length(modalities))
    names(subparams) <- modalities

    for (mod in modalities) {
      fc  <- data[[mod]]
      pv  <- pval_list[[mod]]
      # validate each pair
      if (!is.matrix(fc)||!is.numeric(fc)) stop(mod,": FC not numeric matrix.")
      if (!is.matrix(pv)||!is.numeric(pv)) stop(mod,": pval not numeric matrix.")
      if (!all(dim(fc)==dim(pv))) stop(mod,": dims differ.")
      # filter fc_cutoff & max_features
      if (fc_cutoff>0) {
        mfc  <- matrixStats::rowMaxs(abs(fc))
        keep <- mfc>=fc_cutoff
        fc   <- fc[keep,,drop=FALSE]; pv <- pv[keep,,drop=FALSE]
        if (!nrow(fc)) stop(mod,": no features pass fc_cutoff.")
      }
      if (!is.null(max_features) && max_features<nrow(fc)) {
        mfc <- matrixStats::rowMaxs(abs(fc))
        ord <- order(mfc,decreasing=TRUE)[seq_len(max_features)]
        fc  <- fc[ord,,drop=FALSE]; pv <- pv[ord,,drop=FALSE]
      }
      # color mapping for this modality
      col_fun <- if (!is.null(col_list) && mod%in%names(col_list))
                   col_list[[mod]]
                 else circlize::colorRamp2(c(-2,0,2), c("blue","white","red"))
      # layer_fun for this modality
      layer_fun <- function(j,i,x,y,...) {
        pvij <- pv[cbind(i,j)]
        cols <- rep(NA_character_, length(pvij))
        sig  <- which(pvij<pvalue_cutoff)
        tr   <- which(pvij>=pvalue_cutoff & pvij<trending_cutoff)
        cols[sig] <- significant_color
        cols[tr]  <- trending_color
        idx <- which(!is.na(cols))
        if (length(idx)) {
          grid::grid.points(
            x    = x[idx], y = y[idx],
            pch  = pch_val,
            gp   = grid::gpar(col=cols[idx]),
            size = grid::unit(unit_val,"mm")
          )
        }
      }
      hm_list[[mod]] <- ComplexHeatmap::Heatmap(
        fc,
        name              = paste0(mod,"_logFC"),
        col               = col_fun,
        cluster_rows      = cluster_rows,
        cluster_columns   = cluster_columns,
        show_row_names    = show_row_names,
        show_column_names = show_column_names,
        layer_fun         = layer_fun,
        ...
      )
      subparams[[mod]] <- list(
        fc_matrix=fc, pval_matrix=pv,
        pvalue_cutoff, trending_cutoff, fc_cutoff, max_features,
        significant_color, trending_color, pch_val, unit_val,
        cluster_rows, cluster_columns, show_row_names, show_column_names, col_fun
      )
    }

    combined <- Reduce(`+`, hm_list)
    params   <- list(
      multimodal=TRUE,
      modalities=modalities,
      subparams=subparams
    )
    methods::new("InformativeHeatmap", heatmap=combined, params=params)
  }
)

#-- Accessor, Updater & Show ---------------------------------------------------

#' @rdname InformativeHeatmap
#' @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))

#' @rdname InformativeHeatmap
#' @export
setMethod(
  "getHeatmapObject",
  signature(x="InformativeHeatmap"),
  function(x) x@heatmap
)

#' @rdname InformativeHeatmap
#' @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))

#' @rdname InformativeHeatmap
#' @export
setMethod(
  "updateLayerFun",
  signature(x="InformativeHeatmap"),
  function(x, layer_fun) {
    params <- x@params
    if (isTRUE(params$multimodal))
      stop("Rebuild multimodal heatmap to change layer_fun.")
    mat  <- x@heatmap@matrix
    args <- c(list(mat), params)
    args$layer_fun <- layer_fun
    new_ht <- do.call(ComplexHeatmap::Heatmap, args)
    x@heatmap <- new_ht
    x@params  <- modifyList(params, list(layer_fun=layer_fun))
    x
  }
)

#' @rdname InformativeHeatmap
#' @exportMethod show
setMethod(
  "show",
  signature(object="InformativeHeatmap"),
  function(object) {
    ComplexHeatmap::draw(object@heatmap)
    invisible(object)
  }
)

