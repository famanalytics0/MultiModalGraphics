#’ @title InformativeHeatmap: Create an Enhanced ComplexHeatmap with Significance Layer
#’ @description
#’ This constructor builds an \code{InformativeHeatmap} S4 object from either:
#’ 1. A numeric matrix and (optional) metadata (DE will be computed internally), or
#’ 2. A precomputed fold‐change matrix + matching p‐value matrix (passed via \code{pvalues=}), or
#’ 3. A \code{MultiAssayExperiment} (clusters samples via iClusterPlus, runs limma per cluster, then builds a fold‐change matrix).
#’ 
#’ It applies optional variance filtering, optional feature capping, and then layers significance/trending points on top of a ComplexHeatmap. 
#’ 
#’ @param data Either:
#’   - A numeric matrix (features × samples) **and** \code{pvalues=NULL} (raw expression/count data); DE is computed internally.
#’   - A numeric matrix (features × features) and \code{pvalues=<matrix>} of same dimension (precomputed FC vs. p).
#’   - A \code{MultiAssayExperiment} object with one or more assays.
#’ @param meta Optional \code{data.frame} of sample metadata (only if \code{data} is a raw matrix).
#’ @param assayName Character: only if \code{data} is a \code{MultiAssayExperiment}, the assay to extract.
#’ @param groupColumn Character: grouping column for DE (only used for raw matrix or MAE workflows).
#’ @param sampleType Character or \code{NULL}: faceting column for raw matrix or MAE workflows.
#’ @param timepoint Character or \code{NULL}: faceting column for raw matrix or MAE workflows.
#’ @param dataType Character: \code{"auto"}, \code{"continuous"}, or \code{"count"}. If \code{"auto"}, inferred from \code{data}.
#’ @param var_quantile Numeric in [0,1]: variance filter threshold. Default = \code{0.75}.
#’ @param pvalue_cutoff Numeric: P‐value threshold for “significance” (colored points). Default = \code{0.05}.
#’ @param trending_cutoff Numeric: P‐value threshold for “trending” (intermediate points). Default = \code{0.1}.
#’ @param fc_cutoff Numeric: absolute \code{log2FC} threshold for filtering before heatmap. Default = \code{0.585}.
#’ @param min_samples Integer: minimum samples per group for DE. Default = \code{3}.
#’ @param max_features Integer or \code{NULL}: cap on # of features for heatmap. Default = \code{NULL}.
#’ @param runClustering Logical: whether to run \code{iClusterPlus} on an MAE (only if \code{data} is MAE). Default = \code{TRUE}.
#’ @param K Integer: # of clusters for \code{iClusterPlus}. Default = \code{3}.
#’ @param lambda Numeric: regularization parameter for \code{iClusterPlus}. Default = \code{0.2}.
#’ @param coef Integer: limma contrast index for DE within clusters. Default = \code{2}.
#’ @param pch_val Integer: plotting character for significance dots. Default = \code{16}.
#’ @param unit_val Numeric: dot size in “mm” units. Default = \code{4}.
#’ @param significant_color Character: color for “significant” points (P < pvalue_cutoff). Default = \code{"black"}.
#’ @param trending_color Character: color for “trending” points (pvalue_cutoff ≤ P < trending_cutoff). Default = \code{"yellow"}.
#’ @param heatmap_scale Character: \code{"expression"} or \code{"logFC"}. Default = \code{"logFC"}.
#’ @param ... Additional arguments passed to \code{ComplexHeatmap::Heatmap()}, including a \code{pvalues=} matrix if \code{data} is precomputed FC.
#’
#’ @return An \code{InformativeHeatmap} S4 object with slots:
#’   - \code{@heatmap}: the underlying \code{ComplexHeatmap::Heatmap} object
#’   - \code{@params}: a list of all constructor arguments (for future introspection or updates).
#’ @examples
#’ \dontrun{
#’ # 1) Raw matrix + metadata workflow:
#’ expr_mat <- matrix(rpois(2000, lambda=10), nrow=100, ncol=20)
#’ rownames(expr_mat) <- paste0("Gene", 1:100)
#’ colnames(expr_mat) <- paste0("S", 1:20)
#’ meta_df <- data.frame(
#’   Group = rep(c("G1","G2"), each=10),
#’   SampleType = rep(c("X","Y"), times=10),
#’   stringsAsFactors = FALSE,
#’   row.names = paste0("S", 1:20)
#’ )
#’ ih <- InformativeHeatmap(
#’   data = expr_mat,
#’   meta = meta_df,
#’   groupColumn = "Group",
#’   sampleType  = "SampleType",
#’   dataType    = "count",
#’   pvalue_cutoff = 0.05,
#’   trending_cutoff = 0.1,
#’   fc_cutoff     = 0.585,
#’   max_features  = 50,
#’   col = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
#’   cluster_rows = TRUE,
#’   cluster_columns = TRUE
#’ )
#’ ComplexHeatmap::draw(getHeatmapObject(ih))
#’
#’ # 2) Precomputed FC + p-value matrices:
#’ fc_mat <- matrix(rnorm(200), nrow=20, ncol=10)
#’ pv_mat <- matrix(runif(200), nrow=20, ncol=10)
#’ rownames(fc_mat) <- paste0("Gene", 1:20)
#’ rownames(pv_mat) <- paste0("Gene", 1:20)
#’ colnames(fc_mat) <- colnames(pv_mat) <- paste0("Contrast", 1:10)
#’ ih2 <- InformativeHeatmap(
#’   data = fc_mat,
#’   pvalues = pv_mat,
#’   pvalue_cutoff = 0.05,
#’   trending_cutoff = 0.1,
#’   col = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick"))
#’ )
#’ ComplexHeatmap::draw(getHeatmapObject(ih2))
#’ }
#’ @export
setGeneric("InformativeHeatmap", function(data, ...) standardGeneric("InformativeHeatmap"))

#’ @rdname InformativeHeatmap
setMethod("InformativeHeatmap", "ANY", function(
  data,
  meta               = NULL,
  assayName          = NULL,
  groupColumn        = NULL,
  sampleType         = NULL,
  timepoint          = NULL,
  dataType           = c("auto", "continuous", "count"),
  var_quantile       = 0.75,
  pvalue_cutoff      = 0.05,
  trending_cutoff    = 0.1,
  fc_cutoff          = 0.585,
  min_samples        = 3,
  max_features       = NULL,
  runClustering      = TRUE,
  K                  = 3,
  lambda             = 0.2,
  coef               = 2,
  pch_val            = 16,
  unit_val           = 4,
  significant_color  = "black",
  trending_color     = "yellow",
  heatmap_scale      = c("logFC", "expression"),
  ...
) {
  # Helper: check if precomputed pvalues passed
  extra_args <- list(...)
  has_pvalues <- "pvalues" %in% names(extra_args)
  pvalues_mat <- NULL
  if (has_pvalues) {
    pvalues_mat <- extra_args$pvalues
    extra_args$pvalues <- NULL
    if (!is.matrix(pvalues_mat)) {
      stop("`pvalues` must be a matrix with same dimensions as `data`.")
    }
    if (!all(dim(pvalues_mat) == dim(data))) {
      stop("Dimensions of `pvalues` must match `data` (fold‐change) matrix.")
    }
  }

  # Case A: data is a MultiAssayExperiment
  if (inherits(data, "MultiAssayExperiment")) {
    mae <- data
    if (is.null(assayName)) {
      stop("When `data` is a MultiAssayExperiment, you must supply `assayName`.")
    }
    if (!assayName %in% names(experiments(mae))) {
      stop("Assay '", assayName, "' not found in MAE.")
    }
    # ---- A1) Extract assay matrices for all modalities in MAE ----
    se_list <- lapply(experiments(mae), function(x) {
      mat <- as.matrix(SummarizedExperiment::assay(x))
      mode(mat) <- "numeric"
      t(mat)  # transpose: samples × features
    })
    # Run iClusterPlus if requested
    if (runClustering) {
      if (!requireNamespace("iClusterPlus", quietly = TRUE)) {
        stop("`iClusterPlus` is required for clustering. Install via BiocManager::install('iClusterPlus').")
      }
      dt_list <- stats::setNames(se_list, paste0("dt", seq_along(se_list)))
      fit_ic <- do.call(iClusterPlus, c(dt_list, list(
        type    = rep("gaussian", length(se_list)),
        K       = K,
        lambda  = rep(lambda, length(se_list)),
        maxiter = 20
      )))
      clusters <- factor(fit_ic$clusters)
      design <- stats::model.matrix(~ 0 + clusters)
      colnames(design) <- levels(clusters)
      # Run limma per assay
      limma_results <- lapply(names(se_list), function(name) {
        assay_data <- t(se_list[[name]])  # back to features × samples
        fit_l <- limma::lmFit(assay_data, design)
        fit_l <- limma::eBayes(fit_l)
        tt_l <- limma::topTable(fit_l, coef = coef, number = Inf)
        list(logFC = tt_l$logFC, p_values = tt_l$P.Value, feature = rownames(tt_l))
      })
      names(limma_results) <- names(se_list)
      # Combine logFC and pvalues across assays
      combined_logFC <- do.call(rbind, lapply(limma_results, function(res) res$logFC))
      combined_pvalues <- do.call(rbind, lapply(limma_results, function(res) res$p_values))
      rownames(combined_logFC) <- unlist(lapply(limma_results, function(res) res$feature))
      rownames(combined_pvalues) <- unlist(lapply(limma_results, function(res) res$feature))
      # Transpose for heatmap (features × contrasts)
      fc_mat <- t(combined_logFC)
      pv_mat <- t(combined_pvalues)
      # Now proceed as “precomputed” pathway below:
      data_fc <- fc_mat
      pvalues_fc <- pv_mat
      # Variance filtering & capping
      if (!is.null(var_quantile) && var_quantile > 0 && var_quantile < 1) {
        rv <- matrixStats::rowVars(data_fc, na.rm = TRUE)
        thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
        keep_feats <- rv >= thr
        if (!any(keep_feats)) stop("No features pass variance filter (MAE workflow).")
        data_fc <- data_fc[keep_feats, , drop = FALSE]
        pvalues_fc <- pvalues_fc[keep_feats, , drop = FALSE]
      }
      if (!is.null(max_features) && max_features < nrow(data_fc)) {
        vv <- matrixStats::rowVars(data_fc, na.rm = TRUE)
        top_idx <- order(vv, decreasing = TRUE)[seq_len(max_features)]
        data_fc <- data_fc[top_idx, , drop = FALSE]
        pvalues_fc <- pvalues_fc[top_idx, , drop = FALSE]
      }
      # Construct layer_fun (vectorized)
      layer_fun <- function(j, i, x, y, w, h, fill) {
        ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
        all_inds <- as.vector(ind_mat)
        feat_idx <- as.vector(row(ind_mat))
        contr_idx <- as.vector(col(ind_mat))
        pvals_vec <- pvalues_fc[cbind(feat_idx, contr_idx)]
        col_vec <- ifelse(
          pvals_vec < pvalue_cutoff, significant_color,
          ifelse(pvals_vec < trending_cutoff, trending_color, NA_character_)
        )
        keep <- !is.na(col_vec)
        grid::grid.points(
          x[all_inds][keep],
          y[all_inds][keep],
          pch = pch_val,
          gp  = grid::gpar(col = col_vec[keep]),
          size = grid::unit(unit_val, "mm")
        )
      }
      hm <- ComplexHeatmap::Heatmap(
        data_fc,
        name            = ifelse(heatmap_scale == "logFC", "logFC", "expression"),
        layer_fun       = layer_fun,
        ...
      )
      return(methods::new("InformativeHeatmap", heatmap = hm, params = list(
        data          = "MAE workflow",
        assayName     = assayName,
        groupColumn   = groupColumn,
        sampleType    = sampleType,
        timepoint     = timepoint,
        dataType      = dataType,
        var_quantile  = var_quantile,
        pvalue_cutoff = pvalue_cutoff,
        trending_cutoff= trending_cutoff,
        fc_cutoff     = fc_cutoff,
        min_samples   = min_samples,
        max_features  = max_features,
        runClustering = runClustering,
        K             = K,
        lambda        = lambda,
        coef          = coef,
        pch_val       = pch_val,
        unit_val      = unit_val,
        significant_color = significant_color,
        trending_color    = trending_color,
        heatmap_scale    = heatmap_scale
      )))
    }

  # Case B: data is a raw matrix + metadata
  } else if (is.matrix(data) && is.data.frame(meta)) {
    expr_mat <- data
    meta_df   <- meta
    # Basic checks
    if (!all(colnames(expr_mat) %in% rownames(meta_df))) {
      stop("Column names of `data` matrix must match row names of `meta`.")
    }
    meta2 <- meta_df[colnames(expr_mat), , drop = FALSE]
    # Drop NAs if requested
    if (!is.null(groupColumn) && dropNA) {
      needed_cols <- c(groupColumn, sampleType, timepoint)
      keep <- !Reduce(`|`, lapply(needed_cols, function(col) is.na(meta2[[col]])))
      if (!all(keep)) {
        warning(sum(!keep), " sample(s) removed due to NA in grouping/faceting columns.")
      }
      expr_mat <- expr_mat[, keep, drop = FALSE]
      meta2 <- meta2[keep, , drop = FALSE]
    }
    # Check group sizes
    if (!groupColumn %in% colnames(meta2)) {
      stop("`groupColumn` '", groupColumn, "' not found in metadata.")
    }
    tbl_gc <- table(meta2[[groupColumn]])
    if (any(tbl_gc < min_samples)) {
      stop("At least one level of '", groupColumn, "' has fewer than ", min_samples, " samples.")
    }
    if (ncol(expr_mat) < min_samples * 2) {
      stop("Too few total samples (", ncol(expr_mat), ") after filtering; need at least ", min_samples * 2, ".")
    }
    # Determine dataType
    dataType <- match.arg(dataType)
    if (dataType == "auto") {
      if (all(expr_mat == floor(expr_mat), na.rm = TRUE) && max(expr_mat, na.rm = TRUE) > 30) {
        dataType <- "count"
      } else {
        dataType <- "continuous"
      }
    }
    # Variance filter
    if (!is.null(var_quantile) && var_quantile > 0 && var_quantile < 1) {
      rv <- matrixStats::rowVars(expr_mat, na.rm = TRUE)
      thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
      keep_feats <- rv >= thr
      if (!any(keep_feats)) stop("No features pass variance filter (raw‐matrix workflow).")
      expr_mat <- expr_mat[keep_feats, , drop = FALSE]
    }
    if (!is.null(max_features) && max_features < nrow(expr_mat)) {
      vv <- matrixStats::rowVars(expr_mat, na.rm = TRUE)
      top_idx <- order(vv, decreasing = TRUE)[seq_len(max_features)]
      expr_mat <- expr_mat[top_idx, , drop = FALSE]
    }
    # Construct “cells” (identical to ClearScatterplot_table logic)
    if (!is.null(timepoint)) {
      unique_tps <- unique(meta2[[timepoint]])
      unique_sts <- unique(meta2[[sampleType]])
      cells <- expand.grid(timePoint = unique_tps, SampleType = unique_sts, stringsAsFactors = FALSE)
    } else {
      unique_sts <- unique(meta2[[sampleType]])
      cells <- data.frame(timePoint = NA_character_, SampleType = unique_sts, stringsAsFactors = FALSE)
    }
    # Determine parallelization
    vectorized <- match.arg(vectorized)
    use_par <- FALSE
    if (runClustering && parallel) {
      if (vectorized == "vectorized" || (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM))) {
        use_par <- TRUE
      }
    }
    # Per-cell DE (same as ClearScatterplot_table)
    run_cell <- function(i) {
      tp <- cells$timePoint[i]
      st <- cells$SampleType[i]
      if (!is.null(timepoint)) {
        idx <- which(meta2[[timepoint]] == tp & meta2[[sampleType]] == st)
      } else {
        idx <- which(meta2[[sampleType]] == st)
      }
      if (length(idx) < min_samples || length(unique(meta2[[groupColumn]][idx])) < 2) return(NULL)
      ce <- expr_mat[, idx, drop = FALSE]
      cm <- meta2[idx, , drop = FALSE]
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
      data.frame(
        log2fc     = tt$logFC,
        negLog10p  = -log10(tt$P.Value),
        regulation = ifelse(tt$logFC > 0, "up", "down"),
        SampleType = st,
        timePoint  = tp,
        stringsAsFactors = FALSE,
        row.names  = rownames(tt)
      )
    }
    # Execute per-cell
    if (use_par) {
      lst <- BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
    } else {
      lst <- lapply(seq_len(nrow(cells)), run_cell)
    }
    pd <- do.call(rbind, lst)
    if (is.null(pd) || nrow(pd) == 0) {
      stop("No DE results to plot for ANY cell.")
    }
    # Now pd is a DE table: call precomputed pathway
    data_fc <- t(matrix(pd$log2fc, ncol = length(unique(pd$SampleType))))  # dummy reshape
    # Actually, best to call InformativeHeatmapFromMAT, but for simplicity:
    fc_mat <- matrix(pd$log2fc, 
                     nrow = length(unique(pd$SampleType)), 
                     byrow = TRUE)
    pvals_mat <- matrix(10^(-pd$negLog10p), 
                        nrow = length(unique(pd$SampleType)), 
                        byrow = TRUE)
    # Use InformativeHeatmapFromMAT:
    return(InformativeHeatmapFromMAT(
      logFC_matrix    = fc_mat,
      pvalues_matrix  = pvals_mat,
      pvalue_cutoff   = pvalue_cutoff,
      trending_cutoff = trending_cutoff,
      pch_val         = pch_val,
      unit_val        = unit_val,
      significant_color = significant_color,
      trending_color    = trending_color,
      col             = extra_args$col %||% circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
      cluster_rows    = extra_args$cluster_rows %||% TRUE,
      cluster_columns = extra_args$cluster_columns %||% TRUE,
      show_row_names  = extra_args$show_row_names %||% FALSE,
      show_column_names = extra_args$show_column_names %||% FALSE
    ))
  }

  # Case C: data is a precomputed fold‐change matrix + pvalues
  if (is.matrix(data) && has_pvalues) {
    fc_mat <- data
    pvals_mat <- pvalues_mat
    # Variance filtering & capping
    if (!is.null(var_quantile) && var_quantile > 0 && var_quantile < 1) {
      rv <- matrixStats::rowVars(fc_mat, na.rm = TRUE)
      thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
      keep_feats <- rv >= thr
      if (!any(keep_feats)) stop("No features pass variance filter (precomputed workflow).")
      fc_mat <- fc_mat[keep_feats, , drop = FALSE]
      pvals_mat <- pvals_mat[keep_feats, , drop = FALSE]
    }
    if (!is.null(max_features) && max_features < nrow(fc_mat)) {
      vv <- matrixStats::rowVars(fc_mat, na.rm = TRUE)
      top_idx <- order(vv, decreasing = TRUE)[seq_len(max_features)]
      fc_mat <- fc_mat[top_idx, , drop = FALSE]
      pvals_mat <- pvals_mat[top_idx, , drop = FALSE]
    }
    # Construct vectorized layer_fun
    layer_fun <- function(j, i, x, y, w, h, fill) {
      ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
      all_inds <- as.vector(ind_mat)
      feat_idx <- as.vector(row(ind_mat))
      contr_idx <- as.vector(col(ind_mat))
      pvals_vec <- pvals_mat[cbind(feat_idx, contr_idx)]
      col_vec <- ifelse(
        pvals_vec < pvalue_cutoff, significant_color,
        ifelse(pvals_vec < trending_cutoff, trending_color, NA_character_)
      )
      keep <- !is.na(col_vec)
      grid::grid.points(
        x[all_inds][keep],
        y[all_inds][keep],
        pch = pch_val,
        gp  = grid::gpar(col = col_vec[keep]),
        size = grid::unit(unit_val, "mm")
      )
    }
    hm <- ComplexHeatmap::Heatmap(
      fc_mat,
      name            = ifelse(heatmap_scale == "logFC", "logFC", "foldChange"),
      layer_fun       = layer_fun,
      ...
    )
    return(methods::new("InformativeHeatmap", heatmap = hm, params = list(
      data          = "precomputed",
      pvalues       = pvals_mat,
      pvalue_cutoff = pvalue_cutoff,
      trending_cutoff= trending_cutoff,
      pch_val       = pch_val,
      unit_val      = unit_val,
      significant_color = significant_color,
      trending_color    = trending_color,
      heatmap_scale     = heatmap_scale,
      additional_args   = list(...)
    )))
  }

  # If none of the above matched
  stop("`InformativeHeatmap` could not interpret `data`. It must be either:\n",
       "  • A MultiAssayExperiment (with `assayName`, etc.),\n",
       "  • A raw matrix + metadata (dataType, groupColumn, etc.),\n",
       "  • Or a precomputed fold‐change matrix + `pvalues=` matrix.")
})

#’ @title Update the layer function of an InformativeHeatmap
#’ @param x An \code{InformativeHeatmap} object.
#’ @param layer_fun A function to use instead of the stored \code{layer_fun}.
#’ @return The updated \code{InformativeHeatmap} object.
#’ @importFrom ComplexHeatmap Heatmap
#’ @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))

#’ @rdname updateLayerFun
setMethod("updateLayerFun", "InformativeHeatmap", function(x, layer_fun) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("`ComplexHeatmap` is required to update the layer function. ",
         "Please install via BiocManager::install('ComplexHeatmap').")
  }
  params <- x@params
  params$layer_fun <- layer_fun
  # Retrieve original matrix and any other stored args
  # We assume x@heatmap@matrix holds the heatmap data
  data_mat <- x@heatmap@matrix
  args <- c(list(data_mat), params)
  new_hm <- do.call(ComplexHeatmap::Heatmap, args)
  x@heatmap <- new_hm
  x@params <- params
  return(x)
})

#’ @title Retrieve the ComplexHeatmap object from InformativeHeatmap
#’ @param x An \code{InformativeHeatmap} object.
#’ @return The underlying \code{ComplexHeatmap::Heatmap} object.
#’ @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))

#’ @rdname getHeatmapObject
setMethod("getHeatmapObject", "InformativeHeatmap", function(x) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("`ComplexHeatmap` is required to retrieve the Heatmap object. ",
         "Please install via BiocManager::install('ComplexHeatmap').")
  }
  return(x@heatmap)
})
