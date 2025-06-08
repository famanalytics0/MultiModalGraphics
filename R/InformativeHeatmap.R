#' -----------------------------------------------------------------------------
#' InformativeHeatmap: A Class for Enhanced & Multimodal Heatmaps
#'
#' @slot heatmap A ComplexHeatmap Heatmap or HeatmapList object
#' @slot params  A list of all parameters used to build the heatmap
#' @exportClass InformativeHeatmap
#' @import methods
#' @importFrom ComplexHeatmap Heatmap
setClass(
  "InformativeHeatmap",
  slots = c(
    heatmap = "Heatmap",
    params  = "list"
  )
)

#' -----------------------------------------------------------------------------
#' Generic Constructor for InformativeHeatmap
#'
#' @param data    Either (1) a numeric matrix, (2) a MultiAssayExperiment, or
#'                (3) a named list of logFC matrices.
#' @param ...     Further arguments passed to the appropriate method.
#' @export
setGeneric(
  "InformativeHeatmap",
  function(data, ...) standardGeneric("InformativeHeatmap")
)

#' -----------------------------------------------------------------------------
#' Method: matrix + metadata → InformativeHeatmap
#'
#' @param data               Numeric matrix (features × samples)
#' @param meta               data.frame of sample metadata (rownames = colnames(data))
#'                          Required if runClustering = FALSE
#' @param runClustering      Logical; if TRUE runs iClusterPlus → clusters
#'                          else uses groupColumn
#' @param groupColumn        Column in meta for DE grouping
#' @param continuous         One of "auto","continuous","count"
#' @param var_quantile       Quantile [0–1] for variance filter; NULL = off
#' @param min_features       If var_quantile=NULL, keep top N by variance; NULL = off
#' @param pvalue_cutoff      p < this for “significant”
#' @param trending_cutoff    p in [pvalue_cutoff, trending_cutoff) → “trending”
#' @param fc_cutoff          Minimum |log₂FC| to keep
#' @param max_features       Cap on # of features after DE
#' @param significant_color  Color for significant overlays
#' @param trending_color     Color for trending overlays
#' @param pch_val            Plotting character
#' @param unit_val           Point size in mm
#' @param K                  # clusters for iClusterPlus
#' @param lambda             Regularization(s) for iClusterPlus
#' @param coef               Coefficient index or name in limma
#' @param BPPARAM            BiocParallelParam for DE
#' @param heatmap_data_scale "logFC" or "expression"
#' @param cluster_rows       Logical
#' @param cluster_columns    Logical
#' @param show_row_names     Logical
#' @param show_column_names  Logical
#' @param col                ColorRamp2 function
#' @param ...                Additional args to ComplexHeatmap::Heatmap()
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
    ## — argument matching & sanity checks — ##
    continuous       <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)
    if (!is.matrix(data) || !is.numeric(data)) stop("`data` must be numeric matrix.")
    if (!runClustering) {
      if (!is.data.frame(meta)) stop("`meta` must be a data.frame when runClustering=FALSE.")
      if (!(groupColumn %in% colnames(meta))) stop("groupColumn not found in meta.")
      if (!all(colnames(data) %in% rownames(meta)))
        stop("colnames(data) must match rownames(meta).")
      meta <- meta[colnames(data), , drop=FALSE]
    }
    ## — dataType auto-detect — ##
    dataType <- if (continuous=="auto") {
      if (all(data==floor(data)) && max(data,na.rm=TRUE)>30) "count" else "continuous"
    } else continuous
    ## — variance filtering — ##
    n_feat <- nrow(data)
    if (!is.null(var_quantile)) {
      vr <- matrixStats::rowVars(data, na.rm=TRUE)
      thr <- stats::quantile(vr, probs=var_quantile, na.rm=TRUE)
      keep <- vr>=thr
      if (!is.null(min_features) && sum(keep)<min_features) {
        ord <- order(vr,decreasing=TRUE)
        keep[] <- FALSE; keep[ord[seq_len(min_features)]] <- TRUE
      }
      data <- data[keep, , drop=FALSE]
      if (!nrow(data)) stop("No features after variance filtering.")
    } else if (!is.null(min_features)) {
      vr <- matrixStats::rowVars(data, na.rm=TRUE)
      ord <- order(vr,decreasing=TRUE)
      data <- data[ord[seq_len(min_features)], , drop=FALSE]
    }
    ## — design → clustering or provided groups — ##
    if (runClustering) {
      if (!requireNamespace("iClusterPlus", quietly=TRUE))
        stop("Install 'iClusterPlus' for clustering.")
      dt    <- t(data)
      icl   <- iClusterPlus::iClusterPlus(dt1=dt, type="gaussian",
                                          K=K, lambda=rep(lambda,1), maxiter=20)
      clusters <- factor(icl$clusters)
      design   <- stats::model.matrix(~0+clusters)
      colnames(design) <- levels(clusters)
    } else {
      grp <- meta[[groupColumn]]
      grp <- if (!is.factor(grp)) factor(grp) else grp
      if (nlevels(grp)<2) stop("groupColumn must have ≥2 levels.")
      design <- stats::model.matrix(~0+grp)
      colnames(design) <- levels(grp)
    }
    ## — DE via limma/voom — ##
    fit <- if (dataType=="count") {
      v   <- limma::voom(data, design, plot=FALSE)
      limma::lmFit(v, design)
    } else limma::lmFit(data, design)
    fit <- limma::eBayes(fit)
    ## — pick coefficient — ##
    if (is.character(coef)) {
      ci <- which(colnames(design)==coef)
      if (length(ci)!=1) stop("coef name not in design.")
      coef_idx <- ci
    } else {
      coef_idx <- as.integer(coef)
      if (coef_idx<1||coef_idx>ncol(design)) stop("coef index out of range.")
    }
    tt <- limma::topTable(fit, coef=coef_idx, number=Inf, sort.by="none")
    keep <- !is.na(tt$logFC) & !is.na(tt$P.Value)
    tt <- tt[keep, , drop=FALSE]
    if (!nrow(tt)) stop("No DE results.")
    if (fc_cutoff>0) {
      keep <- abs(tt$logFC)>=fc_cutoff
      tt   <- tt[keep, , drop=FALSE]
      if (!nrow(tt)) stop("No features pass fc_cutoff.")
    }
    if (!is.null(max_features) && max_features<nrow(tt)) {
      ord <- order(abs(tt$logFC),decreasing=TRUE)
      tt  <- tt[ord[seq_len(max_features)], , drop=FALSE]
    }
    ## — prepare heatmap matrix & p-vector — ##
    if (heatmap_data_scale=="logFC") {
      hm      <- matrix(tt$logFC, nrow=1,
                        dimnames=list("logFC", rownames(tt)))
      p_vector <- tt$P.Value
    } else {
      feats   <- rownames(tt)
      hm      <- data[feats, , drop=FALSE]
      p_vector <- tt$P.Value
    }
    ## — vectorized overlay — ##
    layer_fun <- function(j,i,x,y,w,h,fill) {
      pv <- if (heatmap_data_scale=="logFC") p_vector[j] else p_vector[i]
      cols <- ifelse(pv<pvalue_cutoff, significant_color,
                     ifelse(pv<trending_cutoff, trending_color, NA))
      kp   <- which(!is.na(cols))
      if (length(kp))
        grid::grid.points(x[kp], y[kp],
                          pch=pch_val,
                          gp=grid::gpar(col=cols[kp]),
                          size=grid::unit(unit_val,"mm"))
    }
    ## — build heatmap — ##
    ht <- ComplexHeatmap::Heatmap(
      hm, name="Value", col=col,
      cluster_rows=cluster_rows, cluster_columns=cluster_columns,
      show_row_names=show_row_names, show_column_names=show_column_names,
      layer_fun=layer_fun, ...
    )
    ## — return object — ##
    params <- as.list(environment())
    methods::new("InformativeHeatmap", heatmap=ht, params=params)
  }
)

#' -----------------------------------------------------------------------------
#' Method: MultiAssayExperiment → InformativeHeatmap
#'
#' Delegates extraction to the matrix + metadata method.
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "MultiAssayExperiment"),
  function(data, assayName = NULL, ...) {
    assays_list <- MultiAssayExperiment::experiments(data)
    if (is.null(assayName)) {
      if (length(assays_list)!=1L) stop("Please specify assayName.")
      assayName <- names(assays_list)[1]
    }
    if (!(assayName %in% names(assays_list)))
      stop("Assay '", assayName, "' not present.")
    se   <- assays_list[[assayName]]
    expr <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se),
                          stringsAsFactors=FALSE)
    InformativeHeatmap(data=expr, meta=meta, ...)
  }
)

#' -----------------------------------------------------------------------------
#' Method: Multimodal (list of FC matrices + list of P matrices)
#'
#' @export
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
    fc_list <- data
    nmods   <- names(fc_list)
    ht_list <- vector("list", length(fc_list))
    names(ht_list) <- nmods
    for (mod in nmods) {
      fm <- fc_list[[mod]]
      pm <- pval_list[[mod]]
      # reuse table-only workflow via InformativeHeatmap_table
      ht_list[[mod]] <- InformativeHeatmap_table(
        fc_matrix=fm, pval_matrix=pm,
        pvalue_cutoff=pvalue_cutoff,
        trending_cutoff=trending_cutoff,
        fc_cutoff=fc_cutoff,
        max_features=max_features,
        significant_color=significant_color,
        trending_color=trending_color,
        pch_val=pch_val, unit_val=unit_val,
        cluster_rows=cluster_rows,
        cluster_columns=cluster_columns,
        show_row_names=show_row_names,
        show_column_names=show_column_names,
        col=col_list[[mod]] %||%
            circlize::colorRamp2(c(-2,0,2),c("blue","white","red")),
        ...
      )@heatmap
    }
    combined <- Reduce(`+`, ht_list)
    params   <- list(multimodal=TRUE, modalities=nmods)
    methods::new("InformativeHeatmap", heatmap=combined, params=params)
  }
)

#' -----------------------------------------------------------------------------
#' Show method for InformativeHeatmap
#' @export
setMethod(
  "show", "InformativeHeatmap",
  function(object) {
    ht <- object@heatmap
    ComplexHeatmap::draw(ht)
    invisible(object)
  }
)

#' -----------------------------------------------------------------------------
#' Retrieve the underlying Heatmap or HeatmapList
#' @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))
setMethod(
  "getHeatmapObject", "InformativeHeatmap",
  function(x) x@heatmap
)

#' -----------------------------------------------------------------------------
#' Update the layer_fun of an InformativeHeatmap (single-modality only)
#' @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))
setMethod(
  "updateLayerFun", "InformativeHeatmap",
  function(x, layer_fun) {
    params <- x@params
    if (isTRUE(params$multimodal))
      stop("Cannot generic updateLayerFun on multimodal object; reconstruct per modality.")
    hm_mat <- x@heatmap@matrix
    args   <- params
    args$heatmap_data_scale <- NULL      # remove
    args$meta               <- NULL
    args$data               <- hm_mat
    args$layer_fun          <- layer_fun
    new_ht <- do.call(ComplexHeatmap::Heatmap, args)
    x@heatmap <- new_ht
    x@params  <- modifyList(params, list(layer_fun=layer_fun))
    x
  }
)

