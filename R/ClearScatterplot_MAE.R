#’ @title ClearScatterplot_MAE: Construct a Volcano from a MultiAssayExperiment
#’ @description
#’ Perform per‐cell differential expression (DE) on an assay within a \code{MultiAssayExperiment} and build a \code{ClearScatterplot} object. 
#’ Splits samples by \code{sampleType} and \code{timepoint}, runs limma (voom if count data), and combines all DE results into one table. 
#’ Internally calls \code{ClearScatterplot()} on the combined table.
#’
#’ @param mae A \code{MultiAssayExperiment} object containing one or more assays.
#’ @param assayName Character: the name of the assay within \code{mae} to use.
#’ @param groupColumn Character: column name in \code{colData} used as the grouping variable for DE (e.g. \code{"TumorVsNormal"}).
#’ @param sampleType Character or \code{NULL}: column name in \code{colData} used for X‐axis faceting. If \code{NULL}, facets only by \code{timepoint}.
#’ @param timepoint Character or \code{NULL}: column name in \code{colData} used for Y‐axis faceting. If \code{NULL}, no row facets.
#’ @param dataType Character: one of \code{"auto"}, \code{"continuous"}, \code{"count"}. If \code{"auto"}, the code inspects the assay.
#’ @param vectorized Character: one of \code{"auto"}, \code{"perCell"}, \code{"vectorized"}. If \code{"auto"}, parallelize if #cells > \code{bpworkers(BPPARAM)}.
#’ @param parallel Logical: whether to attempt parallel execution. If \code{FALSE}, forces sequential (\code{"perCell"}).
#’ @param BPPARAM A \code{BiocParallelParam} object (e.g., \code{MulticoreParam(workers=2)}).
#’ @param var_quantile Numeric in [0,1]: variance filter threshold. Features with row‐variance below this quantile are dropped. Default = \code{0.75}.
#’ @param pvalue_cutoff Numeric: P‐value threshold for “significance” (DE). Default = \code{0.05}.
#’ @param fc_cutoff Numeric: absolute \code{log2FC} threshold for calling “up”/“down”. Default = \code{0.585}.
#’ @param min_samples Integer: minimum number of samples in each group per cell. Default = \code{3}.
#’ @param max_features Integer or \code{NULL}: cap on # of features to keep after variance filtering. Default = \code{NULL} (no cap).
#’ @param dropNA Logical: whether to drop samples with \code{NA} in any of \code{groupColumn}, \code{sampleType}, or \code{timepoint}. Default = \code{TRUE}.
#’ @param ... Additional arguments passed to \code{ClearScatterplot()} (e.g., custom thresholds).
#’
#’ @return A \code{ClearScatterplot} object containing a combined DE table in \code{@data} and \code{@plot} = \code{NULL}.
#’ @examples
#’ \dontrun{
#’ library(MultiAssayExperiment)
#’ data("miniACC", package = "MultiAssayExperiment")
#’ cs <- ClearScatterplot_MAE(
#’   mae          = miniACC,
#’   assayName    = "RNASeq2GeneNorm",
#’   groupColumn  = "C1A.C1B",
#’   sampleType   = "pathologic_stage",
#’   timepoint    = "MethyLevel",
#’   dataType     = "auto",
#’   vectorized   = "auto",
#’   parallel     = TRUE,
#’   BPPARAM      = BiocParallel::bpparam(),
#’   var_quantile = 0.75,
#’   pvalue_cutoff= 0.05,
#’   fc_cutoff    = 0.585
#’ )
#’ cs <- createPlot(cs)
#’ show(cs)
#’ }
#’ @export
ClearScatterplot_MAE <- function(
  mae,
  assayName,
  groupColumn,
  sampleType       = NULL,
  timepoint        = NULL,
  dataType         = c("auto", "continuous", "count"),
  vectorized       = c("auto", "perCell", "vectorized"),
  parallel         = TRUE,
  BPPARAM          = BiocParallel::bpparam(),
  var_quantile     = 0.75,
  pvalue_cutoff    = 0.05,
  fc_cutoff        = 0.585,
  min_samples      = 3,
  max_features     = NULL,
  dropNA           = TRUE,
  ...
) {
  # ---- 1) Basic checks ----
  if (!inherits(mae, "MultiAssayExperiment")) {
    stop("`mae` must be a MultiAssayExperiment.")
  }
  if (!assayName %in% names(experiments(mae))) {
    stop("Assay '", assayName, "' not found in `mae`.")
  }
  # Extract SummarizedExperiment
  se_full <- experiments(mae)[[assayName]]
  expr_full <- SummarizedExperiment::assay(se_full)
  if (!is.matrix(expr_full)) {
    stop("Assay '", assayName, "' must be a matrix or SummarizedExperiment with numeric assay.")
  }
  # Ensure numeric
  mode(expr_full) <- "numeric"
  # Get metadata from assay-level colData
  meta_assay <- as.data.frame(SummarizedExperiment::colData(se_full), stringsAsFactors = FALSE)
  # If grouping or faceting columns not all in assay colData, merge with top-level colData via sampleMap
  needed_cols <- c(groupColumn, sampleType, timepoint)
  if (any(!needed_cols %in% colnames(meta_assay))) {
    # Merge with top-level
    smap <- sampleMap(mae)
    colmap <- smap[smap$assay == assayName, c("assay", "colname", "primary")]
    topMeta <- as.data.frame(colData(mae), stringsAsFactors = FALSE)
    rownames(topMeta) <- rownames(colData(mae))
    assay_meta2 <- meta_assay
    assay_meta2$colname <- rownames(assay_meta2)
    combined_meta <- dplyr::left_join(
      assay_meta2,
      data.frame(primary = rownames(topMeta), topMeta, stringsAsFactors = FALSE),
      by = c("colname" = "primary")
    )
    meta <- combined_meta
    rownames(meta) <- meta$colname
    meta$colname <- NULL
  } else {
    meta <- meta_assay
  }
  # ---- 2) Align samples ----
  shared <- intersect(colnames(expr_full), rownames(meta))
  if (length(shared) == 0) {
    stop("No overlapping samples between assay and metadata.")
  }
  expr <- expr_full[, shared, drop = FALSE]
  meta <- meta[shared, , drop = FALSE]
  # ---- 3) Drop NAs in grouping/faceting columns ----
  needed_cols2 <- unique(needed_cols[!is.null(needed_cols)])
  if (length(needed_cols2) > 0 && dropNA) {
    keep <- !Reduce(`|`, lapply(needed_cols2, function(col) is.na(meta[[col]])))
    if (!all(keep)) {
      warning(sum(!keep), " sample(s) removed due to NA in grouping/faceting columns.")
    }
    expr <- expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }
  # ---- 4) Check group sizes ----
  if (!groupColumn %in% colnames(meta)) {
    stop("`groupColumn` '", groupColumn, "' not found in metadata.")
  }
  tbl_gc <- table(meta[[groupColumn]])
  if (any(tbl_gc < min_samples)) {
    stop("At least one level of '", groupColumn, "' has fewer than ", min_samples, " samples.")
  }
  if (ncol(expr) < min_samples * 2) {
    stop("Too few total samples (", ncol(expr), ") after filtering; need at least ", min_samples * 2, ".")
  }
  # ---- 5) Variance filtering ----
  if (!is.null(var_quantile) && var_quantile > 0 && var_quantile < 1) {
    rv <- matrixStats::rowVars(expr, na.rm = TRUE)
    thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
    keep_genes <- rv >= thr
    if (!any(keep_genes)) {
      stop("No features pass the variance quantile filter (var_quantile = ", var_quantile, ").")
    }
    expr <- expr[keep_genes, , drop = FALSE]
  }
  # ---- 6) Cap max_features ----
  if (!is.null(max_features) && is.numeric(max_features) && max_features < nrow(expr)) {
    vv <- matrixStats::rowVars(expr, na.rm = TRUE)
    top_idx <- order(vv, decreasing = TRUE)[seq_len(max_features)]
    expr <- expr[top_idx, , drop = FALSE]
  }
  # ---- 7) Determine dataType ----
  dataType <- match.arg(dataType)
  if (dataType == "auto") {
    if (all(expr == floor(expr), na.rm = TRUE) && max(expr, na.rm = TRUE) > 30) {
      dataType <- "count"
    } else {
      dataType <- "continuous"
    }
  }
  # ---- 8) Construct “cells” for faceting ----
  if (!is.null(timepoint)) {
    unique_tps <- unique(meta[[timepoint]])
    unique_sts <- unique(meta[[sampleType]])
    cells <- expand.grid(
      timePoint  = unique_tps,
      SampleType = unique_sts,
      stringsAsFactors = FALSE
    )
  } else {
    unique_sts <- unique(meta[[sampleType]])
    cells <- data.frame(
      timePoint  = NA_character_,
      SampleType = unique_sts,
      stringsAsFactors = FALSE
    )
  }
  # ---- 9) Prepare for parallelization ----
  vectorized <- match.arg(vectorized)
  use_par <- FALSE
  if (parallel) {
    if (vectorized == "vectorized") {
      use_par <- TRUE
    } else if (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)) {
      use_par <- TRUE
    }
  }
  # ---- 10) Define per‐cell DE function ----
  run_cell <- function(i) {
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]
    if (!is.null(timepoint)) {
      idx <- which(meta[[timepoint]] == tp & meta[[sampleType]] == st)
    } else {
      idx <- which(meta[[sampleType]] == st)
    }
    # Skip if too few samples or only one group present
    if (length(idx) < min_samples || length(unique(meta[[groupColumn]][idx])) < 2) return(NULL)
    ce <- expr[, idx, drop = FALSE]
    cm <- meta[idx, , drop = FALSE]
    # If continuous but large values, transform
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) {
      ce <- log2(ce + 1)
    }
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) {
      warning("Too few samples in cell: ", st, "/", tp, ". Skipping DE.")
      return(NULL)
    }
    if (dataType == "count") {
      vfit <- limma::voom(ce, design, plot = FALSE)
      fit <- limma::lmFit(vfit, design)
    } else {
      fit <- limma::lmFit(ce, design)
    }
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = 2, number = Inf)
    if (nrow(tt) == 0) {
      warning("No DE results for cell: ", st, "/", tp)
      return(NULL)
    }
    df <- data.frame(
      log2fc     = tt$logFC,
      negLog10p  = -log10(tt$P.Value),
      regulation = ifelse(tt$logFC > 0, "up", "down"),
      SampleType = st,
      timePoint  = tp,
      stringsAsFactors = FALSE,
      row.names  = rownames(tt)
    )
    return(df)
  }
  # ---- 11) Execute per‐cell DE (parallel or sequential) ----
  if (use_par) {
    lst <- BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lst <- lapply(seq_len(nrow(cells)), run_cell)
  }
  pd <- do.call(rbind, lst)
  if (is.null(pd) || nrow(pd) == 0) {
    stop("No DE results to plot for ANY cell.")
  }
  # ---- 12) Call table‐only constructor ----
  cs_obj <- ClearScatterplot(
    data            = pd,
    highLog2fc      = fc_cutoff,
    lowLog2fc       = -fc_cutoff,
    negLog10pValue  = -log10(pvalue_cutoff),
    ...
  )
  return(cs_obj)
}

#’ @title Show method for ClearScatterplot_MAE
#’ @param object A \code{ClearScatterplot} object returned by \code{ClearScatterplot_MAE}.
#’ @export
setMethod("show", "ClearScatterplot_MAE", function(object) {
  if (is.null(object@plot)) {
    stop("Plot has not been created. Run `createPlot()` on the ClearScatterplot object.")
  }
  print(object@plot)
  invisible(NULL)
})
