# Only suppress check notes on old R versions
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data", "n1", "n2"))

#' S4 Class: ClearScatterplot
#'
#' Encapsulates data and plotting logic for volcano plots with flexible faceting.
#'
#' @slot data   data.frame with plot data (includes `category` factor: "down" / "neutral" / "up")
#' @slot plot   ANY storing the generated ggplot object
#' @import methods
#' @exportClass ClearScatterplot
setClass(
  "ClearScatterplot",
  slots = c(data = "data.frame", plot = "ANY"),
  validity = function(object) {
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ","))
    # Also require `category` to exist (it will be created in the constructor)
    if (!"category" %in% names(object@data)) {
      stop("Slot `data` must contain a factor column `category` (levels: down/neutral/up).")
    }
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
#' @param data           data.frame with required columns: log2fc, negLog10p, regulation, SampleType
#' @param highLog2fc     threshold for up‐regulation
#' @param lowLog2fc      threshold for down‐regulation
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
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ","))

  # Remove rows with NA in key numeric columns
  na_idx <- is.na(data$log2fc) | is.na(data$negLog10p)
  if (any(na_idx)) {
    warning(sum(na_idx), " rows removed due to NA in log2fc or negLog10p")
    data <- data[!na_idx, , drop = FALSE]
  }

  # Compute color_flag as before (numeric -1/0/1)
  data$color_flag <- with(
    data,
    ifelse(
      log2fc > highLog2fc & negLog10p > negLog10pValue,  1,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1, 0)
    )
  )

  # Create a semantic factor column `category` with levels "down", "neutral", "up"
  # -1 → "down", 0 → "neutral", 1 → "up"
  data$category <- factor(
    data$color_flag,
    levels = c(-1L, 0L, 1L),
    labels = c("down", "neutral", "up"),
    ordered = TRUE
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
  df <- data.frame(
    log2fc     = tt$logFC,
    negLog10p  = -log10(tt$P.Value),
    regulation = ifelse(tt$logFC > 0, "up", "down"),
    SampleType = cell$SampleType,
    timePoint  = cell$timePoint,
    stringsAsFactors = FALSE,
    row.names = rownames(tt)
  )
  return(df)
}

#' Constructor: ClearScatterplot from MAE
#'
#' Handles filtering, NA removal, and optional merging of top‐level metadata.
#'
#' @param mae          MultiAssayExperiment object
#' @param assayName    assay name within mae
#' @param groupColumn  grouping column in colData
#' @param sampleType   facet column in colData
#' @param timepoint    optional row facet column in colData
#' @param dataType     one of "auto","continuous","count"
#' @param vectorized   one of "auto","perCell","vectorized"
#' @param parallel     logical; if FALSE, forces sequential
#' @param BPPARAM      BiocParallelParam
#' @param var_quantile variance filter quantile (0–1)
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
  var_quantile= 0.75
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)

  # Check assay existence
  assays <- names(MultiAssayExperiment::experiments(mae))
  if (!assayName %in% assays) stop("Assay '", assayName, "' not found in MAE")

  se   <- MultiAssayExperiment::experiments(mae)[[assayName]]
  expr <- SummarizedExperiment::assay(se)

  # Feature variance filtering
  rv  <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  # Assay-level metadata
  meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)

  # Merge missing metadata columns from top-level MAE if needed
  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  miss <- setdiff(needed, names(meta))
  if (length(miss)) {
    sm  <- MultiAssayExperiment::sampleMap(mae)
    sm  <- sm[sm$assay == assayName, ]
    pm  <- as.data.frame(SummarizedExperiment::colData(mae), stringsAsFactors = FALSE)
    pm2 <- pm[sm$primary, miss, drop = FALSE]
    rownames(pm2) <- sm$colname
    meta <- cbind(meta, pm2[rownames(meta), , drop = FALSE])
  }

  # Match expression columns and metadata rows
  shared <- intersect(colnames(expr), rownames(meta))
  if (length(shared) == 0) stop("No overlapping samples between expression and metadata")
  expr <- expr[, shared, drop = FALSE]
  meta <- meta[shared, , drop = FALSE]

  # Drop samples with NA in required metadata
  keep <- rowSums(is.na(meta[needed])) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  # Ensure each group has ≥ 3 samples and total ≥ 6 samples
  if (any(table(meta[[groupColumn]]) < 3)) stop("Each group needs ≥3 samples")
  if (ncol(expr) < 6) stop("Too few samples after filtering (need ≥6)")

  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Constructor: ClearScatterplot_table
#'
#' Create a ClearScatterplot object from a plain expression matrix and metadata.
#'
#' @param expr         numeric matrix of features (rows) by samples (columns)
#' @param meta         data.frame of sample metadata; row.names must match column names of `expr`
#' @param groupColumn  name of the grouping column in `meta` used for differential testing
#' @param sampleType   name of the facet column in `meta` for plotting
#' @param timepoint    optional name of an additional facet row column in `meta`; if NULL, all data is a single time point
#' @param dataType     one of "auto", "continuous", or "count"; "auto" infers type based on integer status
#' @param vectorized   one of "auto", "perCell", or "vectorized"; "auto" will parallelize if number of cells exceeds workers
#' @param parallel     logical; if FALSE, disables all parallel execution regardless of `vectorized`
#' @param BPPARAM      a BiocParallelParam object for parallel execution
#' @param var_quantile numeric between 0 and 1; quantile threshold for filtering low‐variance features
#'
#' @return A `ClearScatterplot` S4 object ready for plotting via `show()`
#'
#' @importFrom matrixStats rowVars
#' @import BiocParallel
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
  var_quantile= 0.75
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)
  stopifnot(is.matrix(expr), is.data.frame(meta))

  # Filter features by variance
  rv   <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr  <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  # Subset metadata and drop NA rows
  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  meta2  <- meta[, needed, drop = FALSE]
  keep   <- rowSums(is.na(meta2)) == 0
  expr   <- expr[, keep, drop = FALSE]
  meta2  <- meta2[keep, , drop = FALSE]

  # Ensure each group has ≥ 3 samples and at least 1 feature remains
  if (any(table(meta2[[groupColumn]]) < 3)) stop("Each group needs ≥3 samples")
  if (nrow(expr) < 1) stop("No features left after filtering")

  .ClearScatterplot_core(
    expr, meta2, groupColumn, sampleType,
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
  # Drop samples with NA in required columns
  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  keep <- rowSums(is.na(meta[, cols, drop = FALSE])) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  # Auto-detect dataType if needed
  if (dataType == "auto") {
    dataType <- if (all(expr == floor(expr)) && max(expr, na.rm = TRUE) > 30) "count" else "continuous"
  }

  # Build “cells” (timePoint × SampleType combinations)
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

  # Function to run DE for a single cell
  run_cell <- function(i) {
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]
    idx <- if (!is.null(timepoint)) {
      which(meta[[timepoint]] == tp & meta[[sampleType]] == st)
    } else {
      which(meta[[sampleType]] == st)
    }
    if (length(idx) < 2 || length(unique(meta[[groupColumn]][idx])) < 2) return(NULL)
    ce <- expr[, idx, drop = FALSE]
    cm <- meta[idx, , drop = FALSE]
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) {
      ce <- log2(ce + 1)
    }
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) return(NULL)
    .run_DE(ce, cm, groupColumn, design, dataType, list(SampleType = st, timePoint = tp))
  }

  # Decide whether to run in parallel
  use_par <- (vectorized == "vectorized" || (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
  if (use_par) {
    lst <- BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lst <- lapply(seq_len(nrow(cells)), run_cell)
  }

  pd <- do.call(rbind, lst)
  if (!nrow(pd)) stop("No DE results to plot.")
  rownames(pd) <- NULL

  # Pass through the base constructor, which will add `color_flag` + `category`
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
setMethod("createPlot", "ClearScatterplot", function(object,
  color1          = "cornflowerblue",  # up‐regulated color
  color2          = "grey",            # neutral color
  color3          = "indianred",       # down‐regulated color
  xlab            = expression(log2~fold~change),
  ylab            = expression(-log10~p),
  custom_theme    = NULL,
  point_alpha     = 0.5,
  point_size      = 1.75,
  legend_position = "bottom",
  legend_title    = NULL,
  legend_labels   = NULL,
  text_family     = "sans",
  text_size       = 10,
  ...
) {
  df <- object@data

  # Build faceting formula: use timePoint if it exists
  has_tp <- "timePoint" %in% names(df) && !all(is.na(df$timePoint))
  facet <- if (has_tp && length(unique(df$timePoint)) > 1) {
    "timePoint ~ SampleType"
  } else {
    ". ~ SampleType"
  }

  # 1) Base ggplot call, mapping color to semantic `category` factor
  p <- ggplot2::ggplot(df, ggplot2::aes(
          x = log2fc,
          y = negLog10p,
          color = category     # <-- no numeric magic here
        )) +
       ggplot2::geom_jitter(alpha = point_alpha, size = point_size) +
       ggplot2::geom_point(alpha = point_alpha, size = point_size) +
       ggplot2::labs(x = xlab, y = ylab, color = legend_title, ...) +
       ggplot2::facet_grid(stats::as.formula(facet), space = "free") +
       ggplot2::theme_bw() +
       ggplot2::theme(
         text             = ggplot2::element_text(family = text_family, size = text_size),
         panel.grid.major = ggplot2::element_line(color = "grey80"),
         panel.grid.minor = ggplot2::element_blank(),
         strip.background = ggplot2::element_rect(fill = "white", color = "black"),
         strip.text       = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
         axis.title       = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
         axis.text        = ggplot2::element_text(size = text_size, family = text_family),
         legend.position  = legend_position,
         legend.title     = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
         legend.text      = ggplot2::element_text(size = text_size, family = text_family)
       )

  # 2) Semantic color mapping: “up” → color1, “neutral” → color2, “down” → color3
  col_map <- c(
    "up"      = color1,   # category == "up"
    "neutral" = color2,   # category == "neutral"
    "down"    = color3    # category == "down"
  )
  p <- p + ggplot2::scale_color_manual(
    values = col_map,
    breaks = c("down", "neutral", "up"),  # ensures legend order is (down, neutral, up)
    labels = legend_labels
  )

  # 3) Apply any user‐supplied theme override last
  if (!is.null(custom_theme)) {
    p <- p + custom_theme
  }

  # 4) Add count labels in each facet, grouping by `category`
  if (has_tp && length(unique(df$timePoint)) > 1) {
    up <- df[df$category == "up", ] |>
          dplyr::group_by(timePoint, SampleType) |>
          dplyr::tally(name = "n1")
    dn <- df[df$category == "down", ] |>
          dplyr::group_by(timePoint, SampleType) |>
          dplyr::tally(name = "n2")
  } else {
    up <- df[df$category == "up", ] |>
          dplyr::group_by(SampleType) |>
          dplyr::tally(name = "n1")
    dn <- df[df$category == "down", ] |>
          dplyr::group_by(SampleType) |>
          dplyr::tally(name = "n2")
  }
  if (nrow(up)) {
    p <- p + ggplot2::geom_text(
      data    = up,
      ggplot2::aes(label = n1),
      x       = Inf, y = Inf,
      hjust   = 1.1, vjust = 1.1,
      color   = color1,      # “up” labels in color1
      family  = text_family,
      size    = text_size / ggplot2::.pt
    )
  }
  if (nrow(dn)) {
    p <- p + ggplot2::geom_text(
      data    = dn,
      ggplot2::aes(label = n2),
      x       = -Inf, y = Inf,
      hjust   = -0.1, vjust = 1.1,
      color   = color3,      # “down” labels in color3
      family  = text_family,
      size    = text_size / ggplot2::.pt
    )
  }

  # 5) Store and return invisibly
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














