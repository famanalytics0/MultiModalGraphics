#’ @title InformativeHeatmapFromMAT: Build a Heatmap from Precomputed Fold Changes & P‐values
#’ @description
#’ When you already have a fold‐change matrix \code{logFC_matrix} and a matching p‐value matrix \code{pvalue_matrix},
#’ this function wraps them into an \code{InformativeHeatmap} object without any upstream DE. It will apply
#’ variance filtering, optional feature capping, and layer “significance” dots on top of a ComplexHeatmap.
#’ 
#’ @param logFC_matrix Numeric matrix of fold‐changes (features × contrasts).
#’ @param pvalue_matrix Numeric matrix of p‐values (features × contrasts), same dimensions as \code{logFC_matrix}.
#’ @param pvalue_cutoff Numeric: P‐value threshold for “significance”. Default = \code{0.05}.
#’ @param trending_cutoff Numeric: P‐value threshold for “trending”. Default = \code{0.1}.
#’ @param pch_val Integer: plotting character for points. Default = \code{16}.
#’ @param unit_val Numeric: size of points in “mm”. Default = \code{4}.
#’ @param significant_color Character: color for “significant” points. Default = \code{"black"}.
#’ @param trending_color Character: color for “trending” points. Default = \code{"yellow"}.
#’ @param max_features Integer or \code{NULL}: cap on # of features to keep. Default = \code{NULL}.
#’ @param var_quantile Numeric in [0,1]: variance filter threshold. Default = \code{0.75}.
#’ @param col Color mapping (vector or \code{colorRamp2}). Default = \code{NULL} (uses blue-white-red ramp).
#’ @param cluster_rows Logical: whether to cluster rows. Default = \code{TRUE}.
#’ @param cluster_columns Logical: whether to cluster columns. Default = \code{TRUE}.
#’ @param show_row_names Logical: whether to show row names. Default = \code{FALSE}.
#’ @param show_column_names Logical: whether to show column names. Default = \code{FALSE}.
#’ @param ... Additional arguments passed to \code{ComplexHeatmap::Heatmap()}.
#’ 
#’ @return An \code{InformativeHeatmap} object.
#’ @examples
#’ \dontrun{
#’ fc_mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
#’ pv_mat <- matrix(runif(200), nrow = 20, ncol = 10)
#’ rownames(fc_mat) <- rownames(pv_mat) <- paste0("Feature", 1:20)
#’ colnames(fc_mat) <- colnames(pv_mat) <- paste0("Contrast", 1:10)
#’ ih <- InformativeHeatmapFromMAT(
#’   logFC_matrix   = fc_mat,
#’   pvalue_matrix  = pv_mat,
#’   pvalue_cutoff  = 0.05,
#’   trending_cutoff= 0.1,
#’   pch_val        = 16,
#’   unit_val       = 4,
#’   significant_color = "darkred",
#’   trending_color    = "darkorange",
#’   max_features   = 15,
#’   var_quantile   = 0.5,
#’   col            = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#’   cluster_rows   = TRUE,
#’   cluster_columns= TRUE,
#’   show_row_names = TRUE
#’ )
#’ ComplexHeatmap::draw(getHeatmapObject(ih))
#’ }
#’ @export
InformativeHeatmapFromMAT <- function(
  logFC_matrix,
  pvalue_matrix,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  pch_val          = 16,
  unit_val         = 4,
  significant_color= "black",
  trending_color   = "yellow",
  max_features     = NULL,
  var_quantile     = 0.75,
  col              = NULL,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = FALSE,
  show_column_names= FALSE,
  ...
) {
  if (!is.matrix(logFC_matrix) || !is.matrix(pvalue_matrix)) {
    stop("`logFC_matrix` and `pvalue_matrix` must both be matrices.")
  }
  if (!all(dim(logFC_matrix) == dim(pvalue_matrix))) {
    stop("Dimensions of `logFC_matrix` and `pvalue_matrix` must match.")
  }
  data_fc <- logFC_matrix
  pvals_fc <- pvalue_matrix
  # Variance filtering
  if (!is.null(var_quantile) && var_quantile > 0 && var_quantile < 1) {
    rv <- matrixStats::rowVars(data_fc, na.rm = TRUE)
    thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
    keep_feats <- rv >= thr
    if (!any(keep_feats)) {
      stop("No features pass variance filter (clear mat).")
    }
    data_fc <- data_fc[keep_feats, , drop = FALSE]
    pvals_fc <- pvals_fc[keep_feats, , drop = FALSE]
  }
  # Cap features
  if (!is.null(max_features) && max_features < nrow(data_fc)) {
    vv <- matrixStats::rowVars(data_fc, na.rm = TRUE)
    top_idx <- order(vv, decreasing = TRUE)[seq_len(max_features)]
    data_fc <- data_fc[top_idx, , drop = FALSE]
    pvals_fc <- pvals_fc[top_idx, , drop = FALSE]
  }
  # Default color ramp if none provided
  if (is.null(col)) {
    col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  }
  # Build vectorized layer_fun
  layer_fun <- function(j, i, x, y, w, h, fill) {
    ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
    all_inds <- as.vector(ind_mat)
    feat_idx <- as.vector(row(ind_mat))
    contr_idx <- as.vector(col(ind_mat))
    pvals_vec <- pvals_fc[cbind(feat_idx, contr_idx)]
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
    name             = "logFC",
    col              = col,
    cluster_rows     = cluster_rows,
    cluster_columns  = cluster_columns,
    show_row_names   = show_row_names,
    show_column_names= show_column_names,
    layer_fun        = layer_fun,
    ...
  )
  methods::new("InformativeHeatmap", heatmap = hm, params = list(
    logFC_matrix    = logFC_matrix,
    pvalue_matrix   = pvalue_matrix,
    pvalue_cutoff   = pvalue_cutoff,
    trending_cutoff = trending_cutoff,
    pch_val         = pch_val,
    unit_val        = unit_val,
    significant_color = significant_color,
    trending_color    = trending_color,
    var_quantile      = var_quantile,
    max_features      = max_features,
    col               = col,
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_columns,
    show_row_names    = show_row_names,
    show_column_names = show_column_names,
    additional_args   = list(...)
  ))
}
