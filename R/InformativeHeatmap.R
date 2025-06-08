################################################################################
# InformativeHeatmap Class – Multimodal Integration
################################################################################

# Suppress notes about global variables
utils::globalVariables(c("i", "j", "x", "y", "w", "h", "fill"))

#' -----------------------------------------------------------------------------
#' @title InformativeHeatmap: A Class for Enhanced Heatmaps
#'
#' Encapsulates a ComplexHeatmap::Heatmap (or combined HeatmapList) and its params.
#'
#' @slot heatmap Heatmap object (or HeatmapList) from ComplexHeatmap
#' @slot params  List of all parameters used
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
#' @describeIn InformativeHeatmap Generic constructor
#' @param data   Primary input: matrix, MultiAssayExperiment, or list
#' @param ...    Additional args passed to methods
#' @export
setGeneric(
  "InformativeHeatmap",
  function(data, ...) standardGeneric("InformativeHeatmap")
)

#' -----------------------------------------------------------------------------
#' @describeIn InformativeHeatmap Method for matrix + metadata
#' @param meta               Data frame of sample metadata (rownames == colnames(data))
#' @param runClustering      Logical: run iClusterPlus if TRUE
#' @param groupColumn        Column in meta for grouping (if no clustering)
#' @param continuous         One of "auto","continuous","count"
#' @param var_quantile       Quantile [0,1] for variance filter
#' @param min_features       If var_quantile=NULL, keep top N by variance
#' @param pvalue_cutoff      p–value threshold
#' @param trending_cutoff    "Trending" p–value threshold
#' @param fc_cutoff          Minimum |log2FC| for DE
#' @param max_features       Cap on number of features after DE
#' @param significant_color  Color for p < pvalue_cutoff
#' @param trending_color     Color for p in [pvalue, trending)
#' @param pch_val            Plotting character
#' @param unit_val           Point size in mm
#' @param K, lambda, coef    iClusterPlus and limma args
#' @param BPPARAM            BiocParallelParam
#' @param heatmap_data_scale "logFC" or "expression"
#' @param cluster_rows, cluster_columns Logical for clustering
#' @param show_row_names, show_column_names Logical for names
#' @param col                 Color mapping function (circlize::colorRamp2)
#' @param ...                 Passed to ComplexHeatmap::Heatmap()
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "matrix"),
  function(
    data,
    meta               = NULL,
    runClustering      = FALSE,
    groupColumn        = NULL,
    continuous         = c("auto","continuous","count"),
    var_quantile       = NULL,
    min_features       = NULL,
    pvalue_cutoff      = 0.05,
    trending_cutoff    = 0.1,
    fc_cutoff          = 0,
    max_features       = NULL,
    significant_color  = "red",
    trending_color     = "orange",
    pch_val            = 16,
    unit_val           = 2,
    K                  = 3,
    lambda             = 0.2,
    coef               = 2,
    BPPARAM            = BiocParallel::bpparam(),
    heatmap_data_scale = c("logFC","expression"),
    cluster_rows       = TRUE,
    cluster_columns    = TRUE,
    show_row_names     = FALSE,
    show_column_names  = FALSE,
    col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    ...
  ) {
    # (All of your existing matrix‐method code goes here, unchanged,
    #  including variance filtering, DE via limma, building heatmap_mat, vectorized layer_fun,
    #  constructing ht <- ComplexHeatmap::Heatmap(...), and packaging params.)
    #
    # At the end:
    methods::new("InformativeHeatmap", heatmap = ht, params = params_list)
  }
)

#' -----------------------------------------------------------------------------
#' @describeIn InformativeHeatmap Method for MultiAssayExperiment inputs
#' @param data      MultiAssayExperiment
#' @param assayName Which assay to extract
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "MultiAssayExperiment"),
  function(data, assayName = NULL, ...) {
    assays_list <- MultiAssayExperiment::experiments(data)
    if (is.null(assayName)) {
      if (length(assays_list) != 1L)
        stop("Multiple assays in MAE; please specify assayName.")
      assayName <- names(assays_list)[1]
    }
    if (!(assayName %in% names(assays_list)))
      stop("Assay '", assayName, "' not found.")
    se   <- assays_list[[assayName]]
    expr <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
    # Delegate to matrix‐method
    InformativeHeatmap(
      data = expr,
      meta = meta,
      ...
    )
  }
)

#' -----------------------------------------------------------------------------
#' @describeIn InformativeHeatmap Method for list‐based multimodal inputs
#' @param data      Named list of FC matrices
#' @param ...       Named list of p–value matrices passed via pval_list=
#' @export
setMethod(
  "InformativeHeatmap",
  signature(data = "list"),
  function(data, pval_list, pvalue_cutoff = 0.05, trending_cutoff = 0.1,
           fc_cutoff = 0, max_features = NULL, significant_color = "red",
           trending_color = "orange", pch_val = 16, unit_val = 2,
           cluster_rows = TRUE, cluster_columns = TRUE,
           show_row_names = TRUE, show_column_names = TRUE,
           col_list = NULL, ...) {
    # (Your existing loop‐over‐modalities code goes here, unchanged,
    #  including validation, per‐modality layer_fun_mod, collecting ht_mod,
    #  reducing via '+', packaging params_list.)
    #
    methods::new("InformativeHeatmap", heatmap = combined_ht, params = params_list)
  }
)

#' -----------------------------------------------------------------------------
#' @title getHeatmapObject
#' @rdname InformativeHeatmap
#' @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))

#' @export
setMethod(
  "getHeatmapObject",
  signature(x = "InformativeHeatmap"),
  function(x) {
    x@heatmap
  }
)

#' -----------------------------------------------------------------------------
#' @title updateLayerFun
#' @rdname InformativeHeatmap
#' @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))

#' @export
setMethod(
  "updateLayerFun",
  signature(x = "InformativeHeatmap"),
  function(x, layer_fun) {
    params <- x@params
    if (isTRUE(params$multimodal))
      stop("Rebuild multimodal InformativeHeatmap to change layer_fun.")
    # Single modality: rebuild with new layer_fun
    mat    <- x@heatmap@matrix
    ht_args<- c(list(mat), params)
    ht_args$layer_fun <- layer_fun
    x@heatmap <- do.call(ComplexHeatmap::Heatmap, ht_args)
    x@params  <- modifyList(params, list(layer_fun = layer_fun))
    x
  }
)

#' -----------------------------------------------------------------------------
#' @title show
#' @rdname InformativeHeatmap
#' @exportMethod show
setMethod(
  "show",
  signature(object = "InformativeHeatmap"),
  function(object) {
    ht <- getHeatmapObject(object)
    ComplexHeatmap::draw(ht)
    invisible(object)
  }
)

