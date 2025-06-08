#’ -----------------------------------------------------------------------------
#’ @title InformativeHeatmapFromMAE: Build an InformativeHeatmap from a MultiAssayExperiment
#’ @description
#’ Extracts an assay and its sample metadata from a `MultiAssayExperiment`, then
#’ delegates to the matrix+metadata constructor of `InformativeHeatmap()`.
#’
#’ @param data                A `MultiAssayExperiment` object containing one or more assays.
#’ @param assayName           Character or integer; which assay to use. If `NULL` and
#’                            the MAE has exactly one assay, that one is chosen.
#’ @param runClustering       Logical; if `TRUE` runs iClusterPlus on the assay data
#’                            to define clusters. Otherwise, uses `groupColumn`.
#’ @param groupColumn         Character; column name in sample metadata to define groups.
#’ @param continuous          One of `"auto"`, `"continuous"`, or `"count"`. Defaults to `"auto"`.
#’ @param var_quantile        Numeric in `[0,1]`; quantile for variance filtering. `NULL` = no filter.
#’ @param min_features        Integer; if `var_quantile = NULL`, keep top N features by variance.
#’                            `NULL` = no filtering by count.
#’ @param pvalue_cutoff       Numeric `(0,1]`; threshold for significance. Default = `0.05`.
#’ @param trending_cutoff     Numeric `(0,1]`; threshold for “trending.” Default = `0.1`.
#’ @param fc_cutoff           Numeric ≥ 0; minimum |log₂FC| for DE. Default = `0`.
#’ @param max_features        Integer; after DE filtering, keep at most this many features.
#’                            `NULL` = no cap.
#’ @param significant_color   Color for `p < pvalue_cutoff`. Default = `"red"`.
#’ @param trending_color      Color for `p` in `[pvalue_cutoff, trending_cutoff)`. Default = `"orange"`.
#’ @param pch_val             Integer; plotting character. Default = `16`.
#’ @param unit_val            Numeric; point size in mm. Default = `2`.
#’ @param K                   Integer; number of clusters for iClusterPlus. Default = `3`.
#’ @param lambda              Numeric or vector; regularization for iClusterPlus. Default = `0.2`.
#’ @param coef                Integer or character; which coefficient in limma. Default = `2`.
#’ @param BPPARAM             A `BiocParallelParam` object for parallel DE. Default = `bpparam()`.
#’ @param heatmap_data_scale  One of `"logFC"` or `"expression"`; basis for heatmap. Default = `"logFC"`.
#’ @param cluster_rows        Logical; cluster rows? Default = `TRUE`.
#’ @param cluster_columns     Logical; cluster columns? Default = `TRUE`.
#’ @param show_row_names      Logical; show row names? Default = `FALSE`.
#’ @param show_column_names   Logical; show column names? Default = `FALSE`.
#’ @param col                 A color mapping function (e.g. from `circlize::colorRamp2`).  
#’                            Default = `circlize::colorRamp2(c(-2,0,2), c("blue","white","red"))`.
#’ @param ...                 Additional arguments passed to `ComplexHeatmap::Heatmap()`.
#’
#’ @return An S4 object of class `InformativeHeatmap`.
#’ @export
#’ @importFrom MultiAssayExperiment experiments colData sampleMap
#’ @importFrom SummarizedExperiment assay
#’ @import BiocParallel
#’ @import methods
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
    continuous         <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)

    ## Validate class & extract assay
    if (!inherits(data, "MultiAssayExperiment")) {
      stop("`data` must be a MultiAssayExperiment.")
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

    ## Delegate to the matrix+metadata method
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
