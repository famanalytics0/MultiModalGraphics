# Only suppress check notes on old R versions
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".data","n1","n2"))
}

#' S4 Class: ClearScatterplot
#'
#' Encapsulates data and plotting logic for volcano plots with flexible faceting.
#'
#' @slot data   data.frame with plot data
#' @slot plot   ANY storing the generated ggplot object
#' @import methods
#' @exportClass ClearScatterplot
setClass(
  "ClearScatterplot",
  slots = c(data = "data.frame", plot = "ANY"),
  validity = function(object) {
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=","))
    TRUE
  }
)

#' Generic: ClearScatterplot
#'
#' @export
setGeneric(
  "ClearScatterplot",
  function(data, ...) standardGeneric("ClearScatterplot")
)

#' Constructor: ClearScatterplot
#'
#' @param data           data.frame with required columns
#' @param highLog2fc     threshold for up-regulation
#' @param lowLog2fc      threshold for down-regulation
#' @param negLog10pValue threshold for significance
#' @export
ClearScatterplot <- function(
  data,
  highLog2fc     = 0.585,
  lowLog2fc      = -0.585,
  negLog10pValue = 1.301
) {
  stopifnot(
    is.data.frame(data),
    is.numeric(data$log2fc),
    is.numeric(data$negLog10p)
  )
  req <- c("log2fc", "negLog10p", "regulation", "SampleType")
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse=","))

  # remove rows with NA in key numeric columns
  na_idx <- is.na(data$log2fc) | is.na(data$negLog10p)
  if (any(na_idx)) {
    warning(sum(na_idx), " rows removed due to NA in log2fc or negLog10p")
    data <- data[!na_idx, , drop = FALSE]
  }

  data$color_flag <- with(
    data,
    ifelse(
      log2fc > highLog2fc & negLog10p > negLog10pValue, 1,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1, 0)
    )
  )

  methods::new("ClearScatterplot", data = data, plot = NULL)
}

#' Internal: run DE per cell
.run_DE <- function(expr, meta, groupColumn, design, dataType, cell) {
  fit <- if (dataType == "count") {
    v <- limma::voom(expr, design, plot = FALSE)
    limma::lmFit(v, design)
  } else {
    limma::lmFit(expr, design)
  }
  fit <- limma::eBayes(fit)
  tt  <- limma::topTable(fit, coef = 2, number = Inf)
  if (!nrow(tt)) {
    warning(sprintf("No DE results for %s/%s", cell$SampleType, cell$timePoint))
    return(NULL)
  }

  data.frame(
    log2fc     = tt$logFC,
    negLog10p  = -log10(tt$P.Value),
    regulation = ifelse(tt$logFC > 0, 'up', 'down'),
    SampleType = cell$SampleType,
    timePoint  = cell$timePoint,
    stringsAsFactors = FALSE,
    row.names = rownames(tt)
  )
}

#' Constructor: ClearScatterplot from MAE
#'
#' Handles filtering, NA removal, and merging metadata from the MAE when needed.
#'
#' @param mae          MultiAssayExperiment object
#' @param assayName    assay to use (must exist in mae)
#' @param groupColumn  grouping column name in colData
#' @param sampleType   facet column name in colData
#' @param timepoint    optional facet row name in colData
#' @param dataType     one of "auto","continuous","count"
#' @param vectorized   one of "auto","perCell","vectorized"; in "auto" mode, invokes parallel when cells > workers
#' @param parallel     logical; if FALSE, disables parallel even if `vectorized = "auto"`
#' @param BPPARAM      BiocParallelParam for parallel execution
#' @param var_quantile variance filter quantile (0-1)
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom matrixStats rowVars
#' @importFrom stats quantile model.matrix
#' @export
ClearScatterplot_MAE <- function(
  mae,
  assayName,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto","continuous","count"),
  vectorized  = c("auto","perCell","vectorized"),
  parallel    = TRUE,
  BPPARAM     = BiocParallel::bpparam(),
  var_quantile = 0.75
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)

  assays <- names(MultiAssayExperiment::experiments(mae))
  if (!assayName %in% assays) stop(
    "Assay '", assayName, "' not found in MAE. Available: ",
    paste(assays, collapse=", ")
  )

  se_full   <- MultiAssayExperiment::experiments(mae)[[assayName]]
  expr_full <- SummarizedExperiment::assay(se_full)

  rv   <- matrixStats::rowVars(expr_full, na.rm = TRUE)
  thr  <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr_f <- expr_full[rv >= thr, , drop = FALSE]

  meta_f <- as.data.frame(
    SummarizedExperiment::colData(se_full),
    stringsAsFactors = FALSE
  )
  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  missing <- setdiff(needed, names(meta_f))
  if (length(missing)) {
    sm <- MultiAssayExperiment::sampleMap(mae)
    sm <- sm[sm$assay == assayName, ]
    pm <- as.data.frame(
      SummarizedExperiment::colData(mae),
      stringsAsFactors = FALSE
    )
    pm2 <- pm[sm$primary, missing, drop = FALSE]
    rownames(pm2) <- sm$colname
    meta_f <- cbind(meta_f, pm2[rownames(meta_f), , drop = FALSE])
  }

  shared <- intersect(colnames(expr_f), rownames(meta_f))
  n_e <- ncol(expr_f); n_m <- nrow(meta_f); n_s <- length(shared)
  if (n_s == 0) stop(
    sprintf("No overlapping samples: %d cols vs %d rows", n_e, n_m)
  ) else if (n_s == n_e && n_s == n_m) {
    message(sprintf("All %d samples match", n_s))
  } else {
    message(sprintf(
      "%d of %d expr cols and %d meta rows match", n_s, n_e, n_m
    ))
  }
  expr <- expr_f[, shared, drop = FALSE]
  meta <- meta_f[shared, , drop = FALSE]

  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  missing_cols <- setdiff(cols, names(meta))
  if (length(missing_cols)) stop(
    "Metadata missing columns: ", paste(missing_cols, collapse=", ")
  )

  keep <- rowSums(is.na(meta[, cols, drop = FALSE])) == 0
  if (sum(!keep) / length(keep) > 0.1) warning(
    sum(!keep), " of ", length(keep), " samples dropped due to NA"
  )
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  tg <- table(meta[[groupColumn]])
  if (any(tg < 3)) stop(
    "Groups with <3 samples: ",
    paste(names(tg)[tg < 3], tg[tg < 3], collapse=", ")
  )
  if (ncol(expr) < 6) stop("Too few samples after filtering")

  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Constructor: ClearScatterplot from matrix + metadata
#'
#' @importFrom matrixStats rowVars
#' @export
ClearScatterplot_table <- function(
  expr,
  meta,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto","continuous","count"),
  vectorized  = c("auto","perCell","vectorized"),
  parallel    = TRUE,
  BPPARAM     = BiocParallel::bpparam(),
  var_quantile = 0.75
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)
  stopifnot(is.matrix(expr), is.data.frame(meta))

  rv   <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr  <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  md   <- meta[, cols, drop = FALSE]
  keep <- rowSums(is.na(md)) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  tg <- table(meta[[groupColumn]])
  if (any(tg < 3)) stop(
    "Groups with <3 samples: ",
    paste(names(tg)[tg < 3], tg[tg < 3], collapse=", ")
  )
  if (nrow(expr) < 1) stop("No features left after filtering; cannot proceed")

  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Internal core
.ClearScatterplot_core <- function(
  expr, meta, groupColumn, sampleType,
  timepoint, dataType, vectorized, BPPARAM
) {
  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  keep <- rowSums(is.na(meta[, cols, drop = FALSE])) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  if (dataType == "auto") {
    dataType <- if (all(expr == floor(expr)) && max(expr, na.rm = TRUE) > 30) "count" else "continuous"
  }

  cells <- if (!is.null(timepoint)) {
    expand.grid(
      timePoint  = unique(meta[[timepoint]]),
      SampleType = unique(meta[[sampleType]]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      SampleType = unique(meta[[sampleType]]),
      timePoint  = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  run_cell <- function(i) {
    tp <- cells$timePoint[i]; st <- cells$SampleType[i]
    idx <- if (!is.null(timepoint)) {
      which(meta[[timepoint]] == tp & meta[[sampleType]] == st)
    } else {
      which(meta[[sampleType]] == st)
    }
    if (length(idx) < 2 || length(unique(meta[[groupColumn]][idx])) < 2) return(NULL)
    ce <- expr[, idx, drop = FALSE]; cm <- meta[idx, , drop = FALSE]
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) ce <- log2(ce + 1)
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) return(NULL)
    .run_DE(ce, cm, groupColumn, design, dataType, list(SampleType = st, timePoint = tp))
  }

  use_par <- (vectorized == "vectorized" || (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
  lst     <- if (use_par) BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM) else lapply(seq_len(nrow(cells)), run_cell)
  pd      <- do.call(rbind, lst)
  if (!nrow(pd)) stop("No DE results to plot.")
  rownames(pd) <- NULL
  ClearScatterplot(pd)
}

#' Generic: createPlot
#'
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' Method: createPlot
#'
#' @exportMethod createPlot
#' @import ggplot2
#' @importFrom dplyr group_by tally
setMethod("createPlot", "ClearScatterplot", function(object, color1 = "cornflowerblue", color2 = "grey", color3 = "indianred") {
  df <- object@data
  facet <- if (!all(is.na(df$timePoint)) && length(unique(df$timePoint)) > 1) "timePoint ~ SampleType" else ". ~ SampleType"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = negLog10p, color = factor(color_flag))) +
    ggplot2::geom_point(alpha = 0.5, size = 1.75) +
    ggplot2::geom_jitter() +
    ggplot2::labs(x = expression(log2~fold~change), y = expression(-log10~p)) +
    ggplot2::scale_color_manual(values = c(color1, color2, color3)) +
    ggplot2::facet_grid(stats::as.formula(facet), space = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text        = ggplot2::element_text(size = 12, face = "bold"),
      axis.title        = ggplot2::element_text(size = 12, face = "bold"),
      axis.text         = ggplot2::element_text(size = 10),
      legend.position   = "bottom"
    )
  up <- df[df$color_flag == 1, ] |> dplyr::group_by(timePoint, SampleType) |> dplyr::tally(name = "n1")
  dn <- df[df$color_flag == -1, ] |> dplyr::group_by(timePoint, SampleType) |> dplyr::tally(name = "n2")
  if (nrow(up)) p <- p + ggplot2::geom_text(data = up, ggplot2::aes(label = n1), x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, color = color1)
  if (nrow(dn)) p <- p + ggplot2::geom_text(data = dn, ggplot2::aes(label = n2), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, color = color3)
  object@plot <- p
  invisible(object)
})

#' Method: show
#'
#' @rdname ClearScatterplot
#' @exportMethod show
setMethod("show", "ClearScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})













