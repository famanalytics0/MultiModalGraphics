#’ @title ClearScatterplot_table: Construct a Volcano from an Expression Matrix + Metadata
#’ @description
#’ Perform per‐cell differential expression (DE) on a raw expression/count matrix (\code{expr}) and accompanying metadata (\code{meta}), then build a \code{ClearScatterplot} object. 
#’ Splits samples by \code{sampleType} and \code{timepoint}, runs limma (voom if count data), and combines all DE results into one table. Internally calls \code{ClearScatterplot()} on the combined table.
#’
#’ @param expr A numeric matrix (features × samples). Row names = feature IDs, Column names = sample IDs.
#’ @param meta A \code{data.frame} of sample metadata. Row names = sample IDs matching \code{colnames(expr)}.
#’ @param groupColumn Character: column name in \code{meta} used as grouping variable for DE.
#’ @param sampleType Character or \code{NULL}: column name in \code{meta} used for X‐axis faceting. If \code{NULL}, no X facet.
#’ @param timepoint Character or \code{NULL}: column name in \code{meta} used for Y‐axis faceting. If \code{NULL}, no Y facet.
#’ @param dataType Character: one of \code{"auto"}, \code{"continuous"}, \code{"count"}. If \code{"auto"}, the function inspects \code{expr}.
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
#’ @return A \code{ClearScatterplot} object with \code{@data} containing the combined DE table, and \code{@plot} = \code{NULL}.
#’ @examples
#’ \dontrun{
#’ library(curatedPCaData)
#’ pcad <- getPCa("Taylor")
#’ se_taylor <- pcad$Taylor
#’ expr_mat <- SummarizedExperiment::assay(se_taylor, "counts")
#’ meta_df <- as.data.frame(SummarizedExperiment::colData(se_taylor), stringsAsFactors = FALSE)
#’ cs <- ClearScatterplot_table(
#’   expr = expr_mat,
#’   meta = meta_df,
#’   groupColumn  = "DiseaseStatus",
#’   sampleType   = "GleasonScore",
#’   timepoint    = "Race",
#’   dataType     = "auto",
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
ClearScatterplot_table <- function(
  expr,
  meta,
  groupColumn     = "Group",
  sampleType      = "SampleType",
  timepoint       = NULL,
  dataType        = c("auto", "continuous", "count"),
  vectorized      = c("auto", "perCell", "vectorized"),
  parallel        = TRUE,
  BPPARAM         = BiocParallel::bpparam(),
  var_quantile    = 0.75,
  pvalue_cutoff   = 0.05,
  fc_cutoff       = 0.585,
  min_samples     = 3,
  max_features    = NULL,
  dropNA          = TRUE,
  ...
) {
  # ---- 1) Basic checks ----
  if (!is.matrix(expr)) {
    stop("`expr` must be a matrix (features × samples).")
  }
  if (!is.data.frame(meta)) {
    stop("`meta` must be a data.frame with row names matching colnames(expr).")
  }
  if (!all(colnames(expr) %in% rownames(meta))) {
    stop("Column names of `expr` and row names of `meta` must match exactly.")
  }
  # Align metadata to expression
  meta2 <- meta[colnames(expr), , drop = FALSE]
  expr2 <- expr
  # ---- 2) Drop NAs in grouping/faceting columns ----
  needed_cols <- c(groupColumn, sampleType, timepoint)
  if (length(needed_cols) > 0 && dropNA) {
    keep <- !Reduce(`|`, lapply(needed_cols, function(col) is.na(meta2[[col]])))
    if (!all(keep)) {
      warning(sum(!keep), " sample(s) removed due to NA in grouping/faceting columns.")
    }
    expr2 <- expr2[, keep, drop = FALSE]
    meta2 <- meta2[keep, , drop = FALSE]
  }
  # ---- 3) Check group sizes ----
  if (!groupColumn %in% colnames(meta2)) {
    stop("`groupColumn` '", groupColumn, "' not found in metadata.")
  }
  tbl_gc <- table(meta2[[groupColumn]])
  if (any(tbl_gc < min_samples)) {
    stop("At least one level of '", groupColumn, "' has fewer than ", min_samples, " samples.")
  }
  if (ncol(expr2) < min_samples * 2) {
    stop("Too few total samples (", ncol(expr2), ") after filtering; need at least ", min_samples * 2, ".")
  }
  # ---- 4) Variance filtering ----
  if (!is.null(var_quantile) && var_quantile > 0 && var_quantile < 1) {
    rv <- matrixStats::rowVars(expr2, na.rm = TRUE)
    thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
    keep_genes <- rv >= thr
    if (!any(keep_genes)) {
      stop("No features pass the variance quantile filter (var_quantile = ", var_quantile, ").")
    }
    expr2 <- expr2[keep_genes, , drop = FALSE]
  }
  # ---- 5) Cap max_features ----
  if (!is.null(max_features) && max_features < nrow(expr2)) {
    vv <- matrixStats::rowVars(expr2, na.rm = TRUE)
    top_idx <- order(vv, decreasing = TRUE)[seq_len(max_features)]
    expr2 <- expr2[top_idx, , drop = FALSE]
  }
  # ---- 6) Determine dataType ----
  dataType <- match.arg(dataType)
  if (dataType == "auto") {
    if (all(expr2 == floor(expr2), na.rm = TRUE) && max(expr2, na.rm = TRUE) > 30) {
      dataType <- "count"
    } else {
      dataType <- "continuous"
    }
  }
  # ---- 7) Construct “cells” ----
  if (!is.null(timepoint)) {
    unique_tps <- unique(meta2[[timepoint]])
    unique_sts <- unique(meta2[[sampleType]])
    cells <- expand.grid(
      timePoint  = unique_tps,
      SampleType = unique_sts,
      stringsAsFactors = FALSE
    )
  } else {
    unique_sts <- unique(meta2[[sampleType]])
    cells <- data.frame(
      timePoint  = NA_character_,
      SampleType = unique_sts,
      stringsAsFactors = FALSE
    )
  }
  # ---- 8) Prepare for parallelization ----
  vectorized <- match.arg(vectorized)
  use_par <- FALSE
  if (parallel) {
    if (vectorized == "vectorized") {
      use_par <- TRUE
    } else if (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)) {
      use_par <- TRUE
    }
  }
  # ---- 9) Define per‐cell DE function ----
  run_cell <- function(i) {
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]
    if (!is.null(timepoint)) {
      idx <- which(meta2[[timepoint]] == tp & meta2[[sampleType]] == st)
    } else {
      idx <- which(meta2[[sampleType]] == st)
    }
    if (length(idx) < min_samples || length(unique(meta2[[groupColumn]][idx])) < 2) return(NULL)
    ce <- expr2[, idx, drop = FALSE]
    cm <- meta2[idx, , drop = FALSE]
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
  # ---- 10) Execute per‐cell DE ----
  if (use_par) {
    lst <- BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lst <- lapply(seq_len(nrow(cells)), run_cell)
  }
  pd <- do.call(rbind, lst)
  if (is.null(pd) || nrow(pd) == 0) {
    stop("No DE results to plot for ANY cell.")
  }
  # ---- 11) Call table‐only constructor ----
  cs_obj <- ClearScatterplot(
    data            = pd,
    highLog2fc      = fc_cutoff,
    lowLog2fc       = -fc_cutoff,
    negLog10pValue  = -log10(pvalue_cutoff),
    ...
  )
  return(cs_obj)
}

#’ @title Show method for ClearScatterplot_table
#’ @param object A \code{ClearScatterplot} object from \code{ClearScatterplot_table}.
#’ @export
setMethod("show", "ClearScatterplot_table", function(object) {
  if (is.null(object@plot)) {
    stop("Plot has not been created. Run `createPlot()` on the ClearScatterplot object.")
  }
  print(object@plot)
  invisible(NULL)
})
