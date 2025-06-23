#’ @title MultiModalPlot: Unified “One‐Stop” Multimodal Plot Builder
#’ @description
#’ A high‐level wrapper that accepts any combination of:
#’   1. MultiAssayExperiment objects,
#’   2. Lists of expr (matrix) + meta (data.frame),
#’   3. Precomputed DE tables (data.frames with log2fc, negLog10p, regulation, SampleType, optional timePoint).
#’ 
#’ Dispatches each input to the appropriate builder:
#’   - Volcano: ThresholdedScatterplot_MAE / ThresholdedScatterplot_table / ThresholdedScatterplot  
#’   - Heatmap: AnnotatedHeatmapFromMAE / AnnotatedHeatmap  
#’ Returns either a ggplot (volcano) or ComplexHeatmap (heatmap), combining multiple panels with patchwork or HeatmapList.
#’
#’ @param inputs         Named list of inputs (MAE, list(expr,meta), or DE‐table df)
#’ @param assayNames     Named character: assay for each MAE
#’ @param groupColumns   Named character: grouping column for each input
#’ @param sampleTypes    Named character (optional): X‐facet column for each input
#’ @param timepoints     Named character (optional): Y‐facet column for each input
#’ @param dataType       "auto"/"continuous"/"count" or named vector per input
#’ @param var_quantile   Numeric in [0,1]: variance filter threshold (heatmaps & volcano pipelines)
#’ @param pvalue_cutoff  Numeric: p‐value cutoff for significance (volcano/heatmap)
#’ @param trending_cutoff Numeric: p‐value cutoff for trending (heatmap only)
#’ @param fc_cutoff      Numeric: absolute log2FC threshold
#’ @param max_features   Integer or NULL: cap on # features (heatmaps)
#’ @param parallel       Logical: use parallel where supported
#’ @param BPPARAM        BiocParallelParam, default bpparam()
#’ @param panel_type     "volcano" or "heatmap"
#’ @param ...            Passed to underlying constructors / createPlot / Heatmap
#’ @importFrom rlang %||%
#’ @importFrom BiocParallel bpparam
#’ @importFrom patchwork wrap_plots
#’ @importFrom YourPkg AnnotatedHeatmapFromMAE createPlot getHeatmapObject
#’ @export
MultiModalPlot <- function(
  inputs,
  assayNames     = NULL,
  groupColumns   = NULL,
  sampleTypes    = NULL,
  timepoints     = NULL,
  dataType       = "auto",
  var_quantile   = 0.75,
  pvalue_cutoff  = 0.05,
  trending_cutoff= 0.1,
  fc_cutoff      = 0.585,
  max_features   = NULL,
  parallel       = TRUE,
  BPPARAM        = BiocParallel::bpparam(),
  panel_type     = c("volcano", "heatmap"),
  ...
) {
  panel_type <- match.arg(panel_type)
  if (!is.list(inputs) || is.null(names(inputs))) {
    stop("`inputs` must be a named list of MAE objects, expr/meta lists, or DE‐table data.frames.")
  }
  n_mods <- length(inputs)
  panel_list <- vector("list", n_mods)
  names(panel_list) <- names(inputs)

  # helper from rlang
  `%||%` <- function(x, y) if (is.null(x)) y else x

  for (mod in names(inputs)) {
    x  <- inputs[[mod]]
    st <- sampleTypes[[mod]] %||% NULL
    tp <- timepoints[[mod]]  %||% NULL
    dt <- if (length(dataType)==1) dataType else (dataType[[mod]] %||% "auto")

    # (A) MAE input
    if (inherits(x, "MultiAssayExperiment")) {
      if (is.null(assayNames[[mod]])) {
        stop(sprintf("For MAE '%s', supply assayNames[['%s']].", mod, mod))
      }
      if (is.null(groupColumns[[mod]])) {
        stop(sprintf("For MAE '%s', supply groupColumns[['%s']].", mod, mod))
      }
      if (panel_type == "volcano") {
        panel_list[[mod]] <- ThresholdedScatterplot_MAE(
          mae           = x,
          assayName     = assayNames[[mod]],
          groupColumn   = groupColumns[[mod]],
          sampleType    = st,
          timepoint     = tp,
          dataType      = dt,
          vectorized    = if (parallel) "auto" else "perCell",
          parallel      = parallel,
          BPPARAM       = BPPARAM,
          var_quantile  = var_quantile,
          pvalue_cutoff = pvalue_cutoff,
          fc_cutoff     = fc_cutoff,
          max_features  = max_features,
          ...
        )
      } else {
        panel_list[[mod]] <- AnnotatedHeatmapFromMAE(
          mae             = x,
          assayName       = assayNames[[mod]],
          groupColumn     = groupColumns[[mod]],
          sampleType      = st,
          timepoint       = tp,
          dataType        = dt,
          var_quantile    = var_quantile,
          pvalue_cutoff   = pvalue_cutoff,
          trending_cutoff = trending_cutoff,
          fc_cutoff       = fc_cutoff,
          max_features    = max_features,
          runClustering   = parallel,
          BPPARAM         = BPPARAM,
          ...
        )
      }
    }

    # (B) expr/meta list
    else if (is.list(x) && is.matrix(x$expr) && is.data.frame(x$meta)) {
      if (is.null(groupColumns[[mod]])) {
        stop(sprintf("For matrix '%s', supply groupColumns[['%s']].", mod, mod))
      }
      if (panel_type == "volcano") {
        panel_list[[mod]] <- ThresholdedScatterplot_table(
          expr          = x$expr,
          meta          = x$meta,
          groupColumn   = groupColumns[[mod]],
          sampleType    = st,
          timepoint     = tp,
          dataType      = dt,
          vectorized    = if (parallel) "auto" else "perCell",
          parallel      = parallel,
          BPPARAM       = BPPARAM,
          var_quantile  = var_quantile,
          pvalue_cutoff = pvalue_cutoff,
          fc_cutoff     = fc_cutoff,
          max_features  = max_features,
          ...
        )
      } else {
        panel_list[[mod]] <- AnnotatedHeatmap(
          data            = x$expr,
          meta            = x$meta,
          groupColumn     = groupColumns[[mod]],
          sampleType      = st,
          timepoint       = tp,
          dataType        = dt,
          var_quantile    = var_quantile,
          pvalue_cutoff   = pvalue_cutoff,
          trending_cutoff = trending_cutoff,
          fc_cutoff       = fc_cutoff,
          max_features    = max_features,
          runClustering   = parallel,
          BPPARAM         = BPPARAM,
          ...
        )
      }
    }

    # (C) DE‐table data.frame
    else if (is.data.frame(x) &&
             all(c("log2fc","negLog10p","regulation","SampleType") %in% colnames(x))) {
      if (panel_type != "volcano") {
        stop("DE‐tables only support panel_type = 'volcano'.")
      }
      panel_list[[mod]] <- ThresholdedScatterplot(
        data           = x,
        highLog2fc     = fc_cutoff,
        lowLog2fc      = -fc_cutoff,
        negLog10pValue = -log10(pvalue_cutoff),
        ...
      )
    }

    else {
      stop(sprintf("Unrecognized input type for '%s'.", mod))
    }
  }

  # Combine & return
  if (panel_type == "volcano") {
    if (n_mods == 1) {
      cs <- createPlot(panel_list[[1]], ...)
      return(cs@plot)
    }
    plots <- lapply(panel_list, function(cs) createPlot(cs, ...)@plot)
    return(wrap_plots(plots))
  } else {
    hts <- lapply(panel_list, getHeatmapObject)
    if (n_mods == 1) return(hts[[1]])
    return(Reduce(`+`, hts))
  }
}
