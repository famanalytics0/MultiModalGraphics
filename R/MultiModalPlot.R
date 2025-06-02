#’ @title MultiModalPlot: Unified “One‐Stop” Multimodal Plot Builder
#’ @description
#’ A high‐level wrapper that accepts any combination of:
#’   1. \code{MultiAssayExperiment} objects (with \code{assayName} specified),
#’   2. Lists of \code{expr} (matrix) + \code{meta} (data.frame),
#’   3. Precomputed DE tables (data.frames with \code{log2fc}, \code{negLog10p}, \code{regulation}, \code{SampleType}, optional \code{timePoint}).
#’ 
#’ It dispatches each input to the appropriate builder—\code{ClearScatterplot_MAE}, \code{ClearScatterplot_table}, or \code{ClearScatterplot}—when \code{panel_type = "volcano"}. 
#’ Or to \code{InformativeHeatmap} / \code{InformativeHeatmapFromMAT} when \code{panel_type = "heatmap"}. 
#’ Returns a combined plot object: a single \code{ggplot} (volcano) or \code{ComplexHeatmap::HeatmapList} (heatmap).
#’ 
#’ @param inputs A **named** list of inputs. Each element can be:
#’   - A \code{MultiAssayExperiment}.  
#’   - A list with elements \code{expr = <matrix>, meta = <data.frame>}.  
#’   - A DE‐table \code{data.frame} with columns \code{log2fc}, \code{negLog10p}, \code{regulation}, \code{SampleType}, (optional) \code{timePoint}.
#’ @param assayNames Named character: for each MAE in \code{inputs}, which assay to use.
#’ @param groupColumns Named character: for each input, which column to use as grouping variable.
#’ @param sampleTypes Named character (optional): for each input, which column to use for X‐axis faceting.
#’ @param timepoints Named character (optional): for each input, which column to use for Y‐axis faceting.
#’ @param dataType Either a single character \code{"auto"}/\code{"continuous"}/\code{"count"} (recycled) or a named vector mapping each modality to one of those.
#’ @param var_quantile Numeric in [0,1]: variance filter threshold. Default = \code{0.75}.
#’ @param pvalue_cutoff Numeric: P‐value threshold for “significance”. Default = \code{0.05}.
#’ @param trending_cutoff Numeric: P‐value threshold for “trending”. Default = \code{0.1}.
#’ @param fc_cutoff Numeric: absolute \code{log2FC} threshold for significance. Default = \code{0.585}.
#’ @param max_features Integer or \code{NULL}: cap on # of features for heatmaps. Default = \code{NULL}.
#’ @param parallel Logical: whether to attempt parallel execution. Default = \code{TRUE}.
#’ @param BPPARAM A \code{BiocParallelParam} object (e.g., \code{bpparam()}). Default = \code{bpparam()}.
#’ @param panel_type Character: \code{"volcano"} or \code{"heatmap"}. Default = \code{"volcano"}.
#’ @param ... Additional, panel‐specific arguments passed to underlying constructors:
#’   - For volcano: forwarded to \code{createPlot()} (e.g., \code{color_up}, \code{title}, \code{custom_theme}).
#’   - For heatmap: forwarded to \code{ComplexHeatmap::Heatmap()} (e.g., \code{col}, \code{cluster_rows}, \code{show_row_names}).
#’ 
#’ @return If \code{panel_type = "volcano"}:
#’   - If single modality: a single \code{ggplot} object.
#’   - If multiple: a combined \code{ggplot} via \code{patchwork}.
#’ If \code{panel_type = "heatmap"}:
#’   - If single modality: a \code{ComplexHeatmap::Heatmap}.
#’   - If multiple: a combined \code{HeatmapList} (side‐by‐side panels).
#’ 
#’ @examples
#’ \dontrun{
#’ # 1) Single MAE volcano:
#’ data("miniACC", package="MultiAssayExperiment")
#’ volcano_plot <- MultiModalPlot(
#’   inputs        = list(ACC = miniACC),
#’   assayNames    = c(ACC = "RNASeq2GeneNorm"),
#’   groupColumns  = c(ACC = "C1A.C1B"),
#’   sampleTypes   = c(ACC = "pathologic_stage"),
#’   timepoints    = c(ACC = "MethyLevel"),
#’   panel_type    = "volcano"
#’ )
#’ print(volcano_plot)
#’ 
#’ # 2) Two (expr,meta) inputs heatmap:
#’ se_rna   <- curatedPCaData::getPCa("Taylor")$Taylor
#’ expr_rna <- SummarizedExperiment::assay(se_rna, "counts")
#’ meta_rna <- as.data.frame(SummarizedExperiment::colData(se_rna), stringsAsFactors=FALSE)
#’ inputs <- list(RNA = list(expr = expr_rna, meta = meta_rna))
#’ ht <- MultiModalPlot(
#’   inputs        = inputs,
#’   groupColumns  = c(RNA = "DiseaseStatus"),
#’   sampleTypes   = c(RNA = "GleasonScore"),
#’   dataType      = "auto",
#’   panel_type    = "heatmap",
#’   col           = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
#’   cluster_rows  = TRUE,
#’   cluster_columns = TRUE
#’ )
#’ ComplexHeatmap::draw(ht)
#’ }
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
  # 1) Basic input checks
  if (!is.list(inputs) || is.null(names(inputs))) {
    stop("`inputs` must be a named list of MAE objects, (expr,meta) lists, or DE-tables.")
  }
  n_mods <- length(inputs)
  # 2) Storage for each panel object
  panel_list <- vector("list", n_mods)
  names(panel_list) <- names(inputs)
  # 3) Iterate over each modality
  for (mod in names(inputs)) {
    x <- inputs[[mod]]
    # A) x is a MultiAssayExperiment
    if (inherits(x, "MultiAssayExperiment")) {
      if (is.null(assayNames) || is.null(assayNames[[mod]])) {
        stop(sprintf("For MAE modality '%s', please supply `assayNames[['%s']]`.", mod, mod))
      }
      if (is.null(groupColumns) || is.null(groupColumns[[mod]])) {
        stop(sprintf("For MAE modality '%s', please supply `groupColumns[['%s']]`.", mod, mod))
      }
      st <- if (!is.null(sampleTypes)) sampleTypes[[mod]] else NULL
      tp <- if (!is.null(timepoints)) timepoints[[mod]] else NULL
      dt <- if (is.character(dataType) && length(dataType) == 1) dataType else {
        if (!is.null(dataType[[mod]])) dataType[[mod]] else "auto"
      }
      if (panel_type == "volcano") {
        panel_list[[mod]] <- ClearScatterplot_MAE(
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
        panel_list[[mod]] <- InformativeHeatmap(
          data            = x,
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
    # B) x is list(expr, meta)
    else if (is.list(x) && is.matrix(x$expr) && is.data.frame(x$meta)) {
      if (is.null(groupColumns) || is.null(groupColumns[[mod]])) {
        stop(sprintf("For matrix modality '%s', please supply `groupColumns[['%s']]`.", mod, mod))
      }
      st <- if (!is.null(sampleTypes)) sampleTypes[[mod]] else NULL
      tp <- if (!is.null(timepoints)) timepoints[[mod]] else NULL
      dt <- if (is.character(dataType) && length(dataType) == 1) dataType else {
        if (!is.null(dataType[[mod]])) dataType[[mod]] else "auto"
      }
      if (panel_type == "volcano") {
        panel_list[[mod]] <- ClearScatterplot_table(
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
        panel_list[[mod]] <- InformativeHeatmap(
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
    # C) x is a DE‐table data.frame
    else if (is.data.frame(x) &&
             all(c("log2fc", "negLog10p", "regulation", "SampleType") %in% colnames(x))) {
      if (panel_type != "volcano") {
        stop(sprintf("DE‐tables only support `panel_type = 'volcano'`; called with '%s'.", panel_type))
      }
      panel_list[[mod]] <- ClearScatterplot(
        data            = x,
        highLog2fc      = fc_cutoff,
        lowLog2fc       = -fc_cutoff,
        negLog10pValue  = -log10(pvalue_cutoff),
        ...
      )
    } else {
      stop(sprintf(
        "Input '%s' is not recognized. Must be an MAE, a list(expr,meta), or DE‐table data.frame.",
        mod
      ))
    }
  }  # end for(mod)

  # 4) Combine panels
  if (panel_type == "volcano") {
    # If only one modality: return its ggplot directly
    if (n_mods == 1) {
      cs_obj <- panel_list[[1]]
      cs_obj <- createPlot(cs_obj, ...)
      return(cs_obj@plot)
    }
    # Else: combine multiple ggplots via patchwork
    library(patchwork)
    gglist <- lapply(panel_list, function(cs_obj) {
      cs2 <- createPlot(cs_obj, ...)
      cs2@plot
    })
    combined <- Reduce(`+`, gglist)
    return(combined)
  } else {
    # panel_type = "heatmap"
    ht_list <- lapply(panel_list, getHeatmapObject)
    if (n_mods == 1) {
      return(ht_list[[1]])
    }
    combined_ht <- Reduce(`+`, ht_list)
    return(combined_ht)
  }
}
