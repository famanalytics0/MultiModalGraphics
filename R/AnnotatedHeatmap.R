################################################################################
# AnnotatedHeatmap.R
################################################################################

#' @title AnnotatedHeatmap: A Class for Enhanced Heatmaps
#' @description
#' Encapsulates a ComplexHeatmap::Heatmap (or combined HeatmapList) together with
#' all parameters used to build it. Supports:
#'   - single‐modality (expression matrix + metadata)
#'   - table‐based FC/p‐value workflows
#'   - one‐or‐many matrix/metadata workflows (side‐by‐side or aggregated)
#'   - multi‐MAE workflows (list of MultiAssayExperiment)
#'   - multimodal integration (named lists of FC/p‐value matrices)
#'
#' @slot heatmap A ComplexHeatmap::Heatmap or HeatmapList
#' @slot params  A list of all parameters used (including sub‐lists)
#' @exportClass AnnotatedHeatmap
#' @importMethodsFrom methods setClass new
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom matrixStats rowVars rowMaxs
#' @importFrom grid unit gpar grid.points
#' @importFrom stats model.matrix quantile
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom BiocParallel bpparam bpworkers bplapply
setClass(
  "AnnotatedHeatmap",
  slots = c(
    heatmap = "Heatmap",
    params  = "list"
  )
)

#— GENERICS -------------------------------------------------------------------

#' @rdname AnnotatedHeatmap
#' @export
setGeneric("AnnotatedHeatmap",
  function(data, meta = NULL, pval_list = NULL, ...) standardGeneric("AnnotatedHeatmap")
)

#' @rdname AnnotatedHeatmap
#' @export
setGeneric("AnnotatedHeatmap_table", function(fc_matrix, pval_matrix, ...) standardGeneric("AnnotatedHeatmap_table"))

#— SINGLE-MATRIX METHOD -----------------------------------------------------

#' @describeIn AnnotatedHeatmap Build from expression matrix + metadata.
#' @exportMethod AnnotatedHeatmap
setMethod(
  "AnnotatedHeatmap",
  signature(data = "matrix"),
  function(
    data, meta = NULL,
    runClustering = FALSE, groupColumn,
    continuous = c("auto","continuous","count"),
    var_quantile = NULL, min_features = NULL,
    pvalue_cutoff = 0.05, trending_cutoff = 0.1,
    fc_cutoff = 0, max_features = NULL,
    significant_color = "red", trending_color = "orange",
    pch_val = 16, unit_val = 2,
    K = 3, lambda = 0.2, coef = 2,
    BPPARAM = BiocParallel::bpparam(),
    heatmap_data_scale = c("logFC","expression"),
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    col = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    legend_title = if (heatmap_data_scale[1]=="logFC") "logFC" else "Value",
    ...
  ) {
    continuous <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)

    ## 1) Differential expression table
    de_res <- .compute_DE_dataframe(
      expr = data, meta = meta,
      runClustering = runClustering, groupColumn = groupColumn,
      continuous = continuous, var_quantile = var_quantile, min_features = min_features,
      pvalue_cutoff = pvalue_cutoff, trending_cutoff = trending_cutoff,
      fc_cutoff = fc_cutoff, max_features = max_features,
      K = K, lambda = lambda, coef = coef, BPPARAM = BPPARAM
    )

    ## 2) Build matrix + layer_fun
    hm_obj <- .make_heatmap_matrix(de_res$de_df, de_res$expr_filt, heatmap_data_scale)
    layer_fn <- .make_layer_fun(
      pvals = hm_obj$pvals,
      pvalue_cutoff = pvalue_cutoff,
      trending_cutoff = trending_cutoff,
      sig_col = significant_color,
      trend_col = trending_color,
      pch_val = pch_val,
      unit_val = unit_val
    )

    ## 3) Call Heatmap (strip any user‐supplied `name`)
    dots <- list(...); dots$name <- NULL
    ht <- do.call(
      Heatmap,
      c(
        list(
          hm_obj$mat,
          name = legend_title,
          col = col,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          show_row_names = show_row_names,
          show_column_names = show_column_names,
          layer_fun = layer_fn
        ),
        dots
      )
    )

    ## 4) Wrap and return
    params <- as.list(environment())[setdiff(names(environment()), c("..."))]
    new("AnnotatedHeatmap", heatmap = ht, params = params)
  }
)

#— TABLE-BASED METHOD ---------------------------------------------------------

#' @describeIn AnnotatedHeatmap Build from precomputed FC & P-value matrices.
#' @exportMethod AnnotatedHeatmap_table
setMethod(
  "AnnotatedHeatmap_table",
  signature(fc_matrix = "matrix", pval_matrix = "matrix"),
  function(
    fc_matrix, pval_matrix,
    pvalue_cutoff = 0.05, trending_cutoff = 0.1,
    fc_cutoff = 0, max_features = NULL,
    significant_color = "red", trending_color = "orange",
    pch_val = 16, unit_val = 2,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE,
    col = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    legend_title = "logFC",
    ...
  ) {
    ## 1) Filter & cap
    filt <- .filter_fc_pval(fc_matrix, pval_matrix, fc_cutoff, max_features)

    ## 2) Layer function
    layer_tab <- .make_layer_fun_tab(
      pvals = filt$pv,
      pvalue_cutoff = pvalue_cutoff,
      trending_cutoff = trending_cutoff,
      sig_col = significant_color,
      trend_col = trending_color,
      pch_val = pch_val,
      unit_val = unit_val
    )

    ## 3) Draw (strip user name)
    dots <- list(...); dots$name <- NULL
    ht <- do.call(
      Heatmap,
      c(
        list(
          filt$fc,
          name = legend_title,
          col = col,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          show_row_names = show_row_names,
          show_column_names = show_column_names,
          layer_fun = layer_tab
        ),
        dots
      )
    )

    params <- as.list(environment())[setdiff(names(formals()), c("fc_matrix","pval_matrix","..."))]
    new("AnnotatedHeatmap", heatmap = ht, params = params)
  }
)

#— MULTI-/SINGLE-MATRIX METHOD ---------------------------------------------

#' @describeIn AnnotatedHeatmap Build from one-or-more expression+metadata pairs:
#'   - length = 1 → side-by-side expr vs DE heatmaps  
#'   - length >1  → aggregate FC/P and draw combined heatmap
#' @exportMethod AnnotatedHeatmap
setMethod(
  "AnnotatedHeatmap",
  signature(data = "list", meta = "list"),
  function(
    data, meta,
    runClustering = FALSE, groupColumn,
    continuous = c("auto","continuous","count"),
    var_quantile = NULL, min_features = NULL,
    pvalue_cutoff = 0.05, trending_cutoff = 0.1,
    fc_cutoff = 0, max_features = NULL,
    K = 3, lambda = 0.2, coef = 2, BPPARAM = BiocParallel::bpparam(),
    heatmap_data_scale = c("logFC","expression"),
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    col = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    legend_title = if (heatmap_data_scale[1]=="logFC") "logFC" else "Value",
    top_n = NULL, sort_by = c("pvalue","abslogFC"),
    ...
  ) {
    continuous <- match.arg(continuous)
    heatmap_data_scale <- match.arg(heatmap_data_scale)
    sort_by <- match.arg(sort_by)

    if (!identical(sort(names(data)), sort(names(meta))))
      stop("Names of data and meta must match.")

    # Run limma for each matrix
    run_one <- function(expr_i, meta_i) {
      .compute_DE_dataframe(
        expr = expr_i, meta = meta_i,
        runClustering = runClustering, groupColumn = groupColumn,
        continuous = continuous, var_quantile = var_quantile, min_features = min_features,
        pvalue_cutoff = pvalue_cutoff, trending_cutoff = trending_cutoff,
        fc_cutoff = fc_cutoff, max_features = max_features,
        K = K, lambda = lambda, coef = coef, BPPARAM = BPPARAM
      )$de_df
    }

    nmat <- length(data)
    # Single‐matrix special case
    if (nmat == 1L) {
      expr1 <- data[[1L]]
      meta1 <- meta[[1L]]
      de_df <- run_one(expr1, meta1)

      # Subset top_n if requested
      if (!is.null(top_n)) {
        ord_feats <- if (sort_by == "pvalue") {
          order(de_df$P.Value)
        } else {
          order(abs(de_df$logFC), decreasing = TRUE)
        }
        sel <- rownames(de_df)[seq_len(min(top_n, nrow(de_df)))]
        de_df <- de_df[sel, , drop = FALSE]
        expr1 <- expr1[sel, , drop = FALSE]
      } else {
        expr1 <- expr1[rownames(de_df), , drop = FALSE]
      }

      # Build FC and P‐value matrices (features × 1)
      fc_mat <- matrix(de_df$logFC, nrow = nrow(de_df),
                       dimnames = list(rownames(de_df), ""))
      pv_mat <- matrix(de_df$P.Value, nrow = nrow(de_df),
                       dimnames = list(rownames(de_df), ""))

      # Expression heatmap
      ht1 <- Heatmap(
        expr1,
        name = "Expression",
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        ...
      )
      # DE heatmap with overlay
      layer_fn <- .make_layer_fun_tab(
        pvals = pv_mat,
        pvalue_cutoff = pvalue_cutoff,
        trending_cutoff = trending_cutoff,
        sig_col = "red",
        trend_col = "orange",
        pch_val = 16,
        unit_val = 2
      )
      ht2 <- Heatmap(
        fc_mat,
        name = legend_title,
        col = col,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        layer_fun = layer_fn,
        ...
      )
      combo <- ht1 + ht2
      return(new("AnnotatedHeatmap", heatmap = combo, params = list(single = TRUE)))
    }

    # Multi‐matrix aggregation
    de_list <- mapply(run_one, data, meta, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    all_feats <- unique(unlist(lapply(de_list, rownames)))
    # Build big FC/P matrices
    fc_mat <- vapply(de_list, function(df) {
      v <- setNames(rep(NA_real_, length(all_feats)), all_feats)
      v[rownames(df)] <- df$logFC; v
    }, numeric(length(all_feats)))
    pv_mat <- vapply(de_list, function(df) {
      v <- setNames(rep(NA_real_, length(all_feats)), all_feats)
      v[rownames(df)] <- df$P.Value; v
    }, numeric(length(all_feats)))
    colnames(fc_mat) <- colnames(pv_mat) <- names(de_list)
    rownames(fc_mat) <- rownames(pv_mat) <- all_feats

    # Delegate to table‐based
    AnnotatedHeatmap_table(
      fc_matrix = fc_mat,
      pval_matrix = pv_mat,
      pvalue_cutoff = pvalue_cutoff,
      trending_cutoff = trending_cutoff,
      fc_cutoff = fc_cutoff,
      max_features = max_features,
      significant_color = "red",
      trending_color = "orange",
      pch_val = 16,
      unit_val = 2,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      col = col,
      legend_title = legend_title,
      ...
    )
  }
)

#— MULTI-MAE LIST METHOD -----------------------------------------------------

#' @describeIn AnnotatedHeatmap
#' Build from multiple MultiAssayExperiment objects:
#'   Extracts each assay (via `assayNames[[mod]]`) + colData,
#'   then calls the list-of-matrices method.
#' @param data        Named list of MultiAssayExperiment
#' @param assayNames  Named character: which assay in each MAE to use
#' @param groupColumns Named character: grouping column in each colData
#' @exportMethod AnnotatedHeatmap
setMethod(
  "AnnotatedHeatmap",
  signature(data = "list", meta = "missing"),
  function(data, assayNames, groupColumns, ...) {
    # sanity
    if (!all(vapply(data, inherits, logical(1), what = "MultiAssayExperiment")))
      stop("When calling AnnotatedHeatmap(list), all elements must be MultiAssayExperiment.")
    if (is.null(assayNames) || is.null(groupColumns))
      stop("You must supply both `assayNames` and `groupColumns` for a list of MAEs.")
    mods <- names(data)
    if (!all(mods %in% names(assayNames)) || !all(mods %in% names(groupColumns)))
      stop("Names of `assayNames` and `groupColumns` must match names(data).")
    # Convert each MAE → expr + meta
    data_list <- setNames(vector("list", length(mods)), mods)
    meta_list <- setNames(vector("list", length(mods)), mods)
    for (mod in mods) {
      mae <- data[[mod]]; an <- assayNames[[mod]]
      assays <- MultiAssayExperiment::experiments(mae)
      if (!an %in% names(assays))
        stop("Assay '", an, "' not found in MAE '", mod, "'.")
      se <- assays[[an]]
      data_list[[mod]] <- SummarizedExperiment::assay(se)
      meta_list[[mod]] <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
    }
    # Delegate to list-of-matrices
    AnnotatedHeatmap(
      data = data_list,
      meta = meta_list,
      ...,
      groupColumn = groupColumns,
      BPPARAM = BiocParallel::bpparam()
    )
  }
)

#— MULTIMODAL FC/PVAL LIST METHOD -------------------------------------------

#' @describeIn AnnotatedHeatmap Integrate named FC/p‐value matrices.
#' @exportMethod AnnotatedHeatmap
setMethod(
  "AnnotatedHeatmap",
  signature(data = "list", pval_list = "list"),
  function(data, pval_list, ...) {
    if (!identical(sort(names(data)), sort(names(pval_list))))
      stop("Names of `data` and `pval_list` must match.")
    hms <- Map(function(fc, pv) {
      AnnotatedHeatmap_table(
        fc_matrix = fc,
        pval_matrix = pv,
        legend_title = NULL,
        ...
      )
    }, data, pval_list)
    combo <- Reduce(`+`, lapply(hms, slot, "heatmap"))
    new("AnnotatedHeatmap", heatmap = combo, params = list(multimodal = TRUE))
  }
)

#— MAE HELPER ---------------------------------------------------------------

#' @title AnnotatedHeatmapFromMAE
#' @description
#' Extracts assay & metadata from a MultiAssayExperiment and calls
#' the matrix+metadata constructor.
#' @param mae       A `MultiAssayExperiment`
#' @param assayName Which assay to use (auto‐selects if only one)
#' @param ...       Passed to `AnnotatedHeatmap(matrix, meta, ...)`
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay colData
#' @export
AnnotatedHeatmapFromMAE <- function(mae, assayName = NULL, ...) {
  if (!inherits(mae, "MultiAssayExperiment"))
    stop("`mae` must be a MultiAssayExperiment.")
  assays_list <- MultiAssayExperiment::experiments(mae)
  if (is.null(assayName)) {
    if (length(assays_list) != 1L)
      stop("MAE contains multiple assays; please specify `assayName`.")
    assayName <- names(assays_list)[1]
  }
  se <- assays_list[[assayName]]
  expr <- SummarizedExperiment::assay(se)
  meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  AnnotatedHeatmap(expr, meta = meta, ...)
}

#— ACCESSORS ---------------------------------------------------------------

#' @rdname AnnotatedHeatmap
#' @export
setGeneric("getHeatmapObject", function(x) standardGeneric("getHeatmapObject"))
setMethod("getHeatmapObject", "AnnotatedHeatmap", function(x) x@heatmap)

#' @rdname AnnotatedHeatmap
#' @export
setGeneric("updateLayerFun", function(x, layer_fun) standardGeneric("updateLayerFun"))
setMethod("updateLayerFun", "AnnotatedHeatmap", function(x, layer_fun) {
  if (isTRUE(x@params$multimodal))
    stop("To update layer_fun for a multimodal object, reconstruct modality by modality.")
  args <- c(list(x@heatmap@matrix), x@params)
  args$layer_fun <- layer_fun
  x@heatmap <- do.call(Heatmap, args)
  x@params <- modifyList(x@params, list(layer_fun = layer_fun))
  x
})

#' @rdname AnnotatedHeatmap
#' @exportMethod show
setMethod("show", "AnnotatedHeatmap", function(object) {
  ComplexHeatmap::draw(object@heatmap)
  invisible(object)
})

#— INTERNAL HELPERS ---------------------------------------------------------

#' @keywords internal
#' @importFrom stats quantile model.matrix
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom matrixStats rowVars
.compute_DE_dataframe <- function(
  expr, meta,
  runClustering, groupColumn,
  continuous, var_quantile, min_features,
  pvalue_cutoff, trending_cutoff,
  fc_cutoff, max_features,
  K, lambda, coef,
  BPPARAM = BiocParallel::bpparam()
) {
  # variance filter
  if (!is.null(var_quantile)) {
    v <- rowVars(expr, na.rm = TRUE)
    thr <- quantile(v, var_quantile, na.rm = TRUE)
    keep <- v >= thr
    if (!is.null(min_features) && sum(keep) < min_features) {
      ord <- order(v, decreasing = TRUE)
      keep[] <- FALSE; keep[ord[seq_len(min_features)]] <- TRUE
    }
    expr <- expr[keep, , drop = FALSE]
  } else if (!is.null(min_features)) {
    v <- rowVars(expr, na.rm = TRUE)
    ord <- order(v, decreasing = TRUE)
    expr <- expr[ord[seq_len(min_features)], , drop = FALSE]
  }

  # design
  if (runClustering) {
    ic <- iClusterPlus::iClusterPlus(
      t(expr), type = "gaussian",
      K = K, lambda = rep(lambda, length = length(lambda)),
      maxiter = 20
    )
    grp <- factor(ic$clusters)
  } else {
    grp <- meta[[groupColumn]]
    if (!is.factor(grp)) grp <- factor(grp)
  }
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)

  # fit
  if (continuous == "count") {
    vfit <- voom(expr, design, plot = FALSE)
    fit <- lmFit(vfit, design)
  } else {
    fit <- lmFit(expr, design)
  }
  fit <- eBayes(fit)
  ci <- if (is.character(coef)) match(coef, colnames(design)) else as.integer(coef)
  tt <- topTable(fit, coef = ci, number = Inf, sort.by = "none")
  tt <- tt[complete.cases(tt[, c("logFC", "P.Value")]), , drop = FALSE]

  # FC & cap
  if (fc_cutoff > 0) tt <- tt[abs(tt$logFC) >= fc_cutoff, , drop = FALSE]
  if (!is.null(max_features) && nrow(tt) > max_features) {
    ord <- order(abs(tt$logFC), decreasing = TRUE)
    tt <- tt[ord[seq_len(max_features)], , drop = FALSE]
  }

  list(expr_filt = expr, de_df = tt)
}

#' @keywords internal
.make_heatmap_matrix <- function(de_df, expr_filt, heatmap_data_scale) {
  if (heatmap_data_scale == "logFC") {
    mat <- matrix(de_df$logFC, nrow = 1, dimnames = list("logFC", rownames(de_df)))
    pvals <- de_df$P.Value
  } else {
    feats <- rownames(de_df)
    mat <- expr_filt[feats, , drop = FALSE]
    pvals <- de_df$P.Value
  }
  list(mat = mat, pvals = pvals)
}

#' @keywords internal
.make_layer_fun <- function(
  pvals, pvalue_cutoff, trending_cutoff,
  sig_col, trend_col, pch_val, unit_val
) {
  force(pvals); force(pvalue_cutoff); force(trending_cutoff)
  force(sig_col); force(trend_col)
  function(j, i, x, y, w, h, fill) {
    pv <- if (length(dim(x)) == 0) pvals[j] else pvals[i]
    cols <- ifelse(pv < pvalue_cutoff, sig_col,
                   ifelse(pv < trending_cutoff, trend_col, NA))
    idx <- which(!is.na(cols))
    if (length(idx)) {
      grid::grid.points(
        x = x[idx], y = y[idx],
        pch = pch_val,
        gp = grid::gpar(col = cols[idx]),
        size = grid::unit(unit_val, "mm")
      )
    }
  }
}

#' @keywords internal
.filter_fc_pval <- function(fc, pv, fc_cutoff, max_features) {
  if (fc_cutoff > 0) {
    keep <- rowMaxs(abs(fc)) >= fc_cutoff
    fc <- fc[keep, , drop = FALSE]
    pv <- pv[keep, , drop = FALSE]
  }
  if (!is.null(max_features) && nrow(fc) > max_features) {
    ord <- order(rowMaxs(abs(fc)), decreasing = TRUE)[seq_len(max_features)]
    fc <- fc[ord, , drop = FALSE]
    pv <- pv[ord, , drop = FALSE]
  }
  list(fc = fc, pv = pv)
}

#' @keywords internal
.make_layer_fun_tab <- function(
  pvals, pvalue_cutoff, trending_cutoff,
  sig_col, trend_col, pch_val, unit_val
) {
  force(pvals); force(pvalue_cutoff); force(trending_cutoff)
  force(sig_col); force(trend_col)
  function(j, i, x, y, w, h, fill) {
    pv <- pvals[cbind(i, j)]
    cols <- rep(NA_character_, length(pv))
    cols[pv < pvalue_cutoff]                         <- sig_col
    cols[pv >= pvalue_cutoff & pv < trending_cutoff]  <- trend_col
    idx <- which(!is.na(cols))
    if (length(idx)) {
      grid::grid.points(
        x = x[idx], y = y[idx],
        pch = pch_val,
        gp = grid::gpar(col = cols[idx]),
        size = grid::unit(unit_val, "mm")
      )
    }
  }
}
