# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2"))

#' ClearScatterplot: S4 Class for Faceted Volcano Plots
#'
#' This S4 class encapsulates the data and plotting logic for volcano (MA) plots,
#' with flexible faceting by sample type and optional timepoint, including overlayed
#' counts of up-/down-regulated features.
#'
#' @slot data       data.frame with columns: log2fc, negLog10p, regulation,
#'                  SampleType, timePoint (optional), color_flag
#' @slot plot       ggplot object (or NULL until created)
#'
#' @details
#' `ClearScatterplot` objects are typically constructed via
#' [ClearScatterplot_MAE()] or [ClearScatterplot_table()], which handle normalization,
#' DE, and input validation. Plots are generated with [createPlot()].
#'
#' @import methods
#' @exportClass ClearScatterplot
setClass(
  "ClearScatterplot",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  ),
  validity = function(object) {
    req_cols <- c("log2fc", "negLog10p", "regulation", "SampleType")
    missing  <- setdiff(req_cols, names(object@data))
    if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
    TRUE
  }
)

#' ClearScatterplot Constructor
#'
#' Validates and builds a ClearScatterplot S4 object from a results data.frame.
#'
#' @param data           data.frame with required columns
#' @param highLog2fc     Threshold for up-regulation (default 0.585)
#' @param lowLog2fc      Threshold for down-regulation (default -0.585)
#' @param negLog10pValue Threshold for significance (-log10 p, default 1.301)
#'
#' @details
#' The `data` frame must contain: log2fc, negLog10p, regulation, SampleType.
#' color_flag (for plot coloring) will be computed automatically.
#'
#' @return S4 ClearScatterplot object
#'
#' @examples
#' df <- data.frame(
#'   log2fc = rnorm(100),
#'   negLog10p = runif(100, 0, 5),
#'   regulation = sample(c("up", "down"), 100, TRUE),
#'   SampleType = sample(c("A", "B"), 100, TRUE),
#'   timePoint = sample(c("T1", "T2"), 100, TRUE)
#' )
#' cs <- ClearScatterplot(df)
#' cs <- createPlot(cs)
#' show(cs)
#' @export
ClearScatterplot <- function(
  data,
  highLog2fc     = 0.585,
  lowLog2fc      = -0.585,
  negLog10pValue = 1.301
) {
  stopifnot(is.data.frame(data))
  req_cols <- c("log2fc", "negLog10p", "regulation", "SampleType")
  missing  <- setdiff(req_cols, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  data$color_flag <- with(data,
    ifelse(log2fc > highLog2fc & negLog10p > negLog10pValue, 1,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1, 0)
    )
  )
  new("ClearScatterplot", data = data, plot = NULL)
}

#' Internal: Differential Expression for a cell (SampleType x TimePoint)
#'
#' Used internally by ClearScatterplot_MAE and ClearScatterplot_table.
#'
#' @noRd
.run_DE <- function(expr, meta, groupColumn, designMat, dataType, cellLabel) {
  if (dataType == "count") {
    v <- limma::voom(expr, designMat, plot = FALSE)
    fit <- limma::lmFit(v, designMat)
  } else {
    fit <- limma::lmFit(expr, designMat)
  }
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf)
  if (nrow(tt) == 0) return(NULL)
  data.frame(
    log2fc     = tt[["logFC"]],
    negLog10p  = -log10(tt[["P.Value"]]),
    regulation = ifelse(tt[["logFC"]] > 0, "up", "down"),
    SampleType = cellLabel$SampleType,
    timePoint  = cellLabel$timePoint,
    stringsAsFactors = FALSE,
    row.names = rownames(tt)
  )
}

#' Build ClearScatterplot from MultiAssayExperiment
#'
#' Runs cell-wise (SampleType x timePoint) differential expression and returns
#' a ClearScatterplot object. Handles RNA-seq (count), microarray (continuous),
#' and automatically detects expression type by default.
#'
#' @param mae         MultiAssayExperiment containing assay and metadata
#' @param assayName   Name of the SummarizedExperiment in experiments(mae)
#' @param groupColumn Column in colData defining the group (e.g., case/control)
#' @param sampleType  Column for faceting columns
#' @param timepoint   Column for faceting rows (optional; NULL for cross-sectional)
#' @param dataType    "auto" (default), "continuous" or "count"
#' @param vectorized  "auto", "perCell" (default per-cell DE), or "vectorized" (full model)
#' @param BPPARAM     BiocParallelParam object for parallelization (default: bpparam())
#'
#' @details
#' - For RNAseq, use "count" (applies voom)
#' - For log-scale or microarray, use "continuous"
#' - "auto" inspects data and picks appropriately
#' - Handles NA and group imbalance gracefully
#'
#' @return S4 ClearScatterplot object
#'
#' @examples
#' library(MultiAssayExperiment)
#' data("miniACC", package = "MultiAssayExperiment")
#' cs <- ClearScatterplot_MAE(
#'   mae         = miniACC,
#'   assayName   = "RNASeq2GeneNorm",
#'   groupColumn = "C1A.C1B",
#'   sampleType  = "pathologic_stage",
#'   timepoint   = "MethyLevel",
#'   dataType    = "auto",
#'   vectorized  = "auto"
#' )
#' cs <- createPlot(cs)
#' show(cs)
#' @export
ClearScatterplot_MAE <- function(
  mae,
  assayName,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto", "continuous", "count"),
  vectorized  = c("auto", "perCell", "vectorized"),
  BPPARAM     = BiocParallel::bpparam()
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)

  if (!inherits(mae, "MultiAssayExperiment")) stop("mae must be a MultiAssayExperiment.")
  assays <- names(MultiAssayExperiment::experiments(mae))
  if (is.null(assayName) || !(assayName %in% assays))
    stop("Specify a valid assayName from: ", paste(assays, collapse = ", "))
  se <- MultiAssayExperiment::experiments(mae)[[assayName]]
  expr <- SummarizedExperiment::assay(se)
  meta <- as.data.frame(MultiAssayExperiment::colData(mae), stringsAsFactors = FALSE)

  .ClearScatterplot_core(
    expr        = expr,
    meta        = meta,
    groupColumn = groupColumn,
    sampleType  = sampleType,
    timepoint   = timepoint,
    dataType    = dataType,
    vectorized  = vectorized,
    BPPARAM     = BPPARAM
  )
}

#' Build ClearScatterplot from Expression Matrix and Metadata Table
#'
#' Same as [ClearScatterplot_MAE()], but works directly with a normalized
#' expression matrix (microarray, voom-counts) and corresponding metadata.
#'
#' @param expr         Gene-by-sample matrix
#' @param meta         data.frame with sample annotations
#' @inheritParams ClearScatterplot_MAE
#'
#' @examples
#' expr <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' meta <- data.frame(
#'   Group = rep(c("A", "B"), 5),
#'   SampleType = rep(c("Tissue", "Blood"), each = 5),
#'   timePoint = rep(c("Day1", "Day2"), each = 5),
#'   row.names = paste0("S", 1:10)
#' )
#' cs <- ClearScatterplot_table(
#'   expr        = expr,
#'   meta        = meta,
#'   groupColumn = "Group",
#'   sampleType  = "SampleType",
#'   timepoint   = "timePoint",
#'   dataType    = "auto"
#' )
#' cs <- createPlot(cs)
#' show(cs)
#' @export
ClearScatterplot_table <- function(
  expr,
  meta,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto", "continuous", "count"),
  vectorized  = c("auto", "perCell", "vectorized"),
  BPPARAM     = BiocParallel::bpparam()
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)
  stopifnot(is.matrix(expr), is.data.frame(meta))
  .ClearScatterplot_core(
    expr        = expr,
    meta        = meta,
    groupColumn = groupColumn,
    sampleType  = sampleType,
    timepoint   = timepoint,
    dataType    = dataType,
    vectorized  = vectorized,
    BPPARAM     = BPPARAM
  )
}

#' @noRd
.ClearScatterplot_core <- function(
  expr, meta, groupColumn, sampleType, timepoint,
  dataType, vectorized, BPPARAM
) {
  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  missing <- setdiff(needed, names(meta))
  if (length(missing)) stop("Missing required metadata: ", paste(missing, collapse = ", "))
  keep <- rownames(meta)[stats::complete.cases(meta[, needed, drop = FALSE])]
  if (length(keep) < 2) stop("Too few samples after NA removal.")
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  if (dataType == "auto") {
    if (all(expr == floor(expr)) && max(expr, na.rm = TRUE) > 30) {
      dataType <- "count"
    } else {
      dataType <- "continuous"
    }
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
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]
    idx <- if (!is.null(timepoint)) {
      which(meta[[timepoint]] == tp & meta[[sampleType]] == st)
    } else {
      which(meta[[sampleType]] == st)
    }
    if (length(idx) < 2 || length(unique(meta[[groupColumn]][idx])) < 2) return(NULL)
    cell_expr <- expr[, idx, drop = FALSE]
    cell_meta <- meta[idx, , drop = FALSE]
    if (dataType == "continuous" && max(cell_expr, na.rm = TRUE) > 50) {
      cell_expr <- log2(cell_expr + 1)
    }
    design <- stats::model.matrix(~ cell_meta[[groupColumn]])
    if (nrow(cell_meta) <= ncol(design)) return(NULL)
    .run_DE(cell_expr, cell_meta, groupColumn, design, dataType, list(
      SampleType = st,
      timePoint  = tp
    ))
  }
  use_parallel <- (vectorized == "vectorized" || (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
  df_list <- if (use_parallel) {
    BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lapply(seq_len(nrow(cells)), run_cell)
  }
  plotdata <- do.call(rbind, df_list)
  if (is.null(plotdata) || nrow(plotdata) == 0) stop("No DE results to plot.")
  rownames(plotdata) <- NULL
  ClearScatterplot(plotdata)
}

#' Render a ClearScatterplot as a ggplot
#'
#' Builds and stores the faceted volcano plot, using theme_bw,
#' facet_grid, and overlays up/down counts per cell.
#'
#' @param object   ClearScatterplot object
#' @param color1   up-regulation color (default "cornflowerblue")
#' @param color2   neutral color (default "grey")
#' @param color3   down-regulation color (default "indianred")
#'
#' @return The input object, with the plot slot populated
#'
#' @examples
#' cs <- ClearScatterplot(df)
#' cs <- createPlot(cs)
#' show(cs)
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))
setMethod(
  "createPlot",
  signature(object = "ClearScatterplot"),
  function(object,
           color1 = "cornflowerblue",
           color2 = "grey",
           color3 = "indianred") {
    df <- object@data
    xvar <- "SampleType"; yvar <- "timePoint"
    ux   <- unique(df[[xvar]]); uy <- unique(df[[yvar]])
    facet_formula <- if (!all(is.na(uy)) && length(uy) > 1 && length(ux) > 1) {
      stats::as.formula(paste(yvar, "~", xvar))
    } else if (!all(is.na(uy)) && length(uy) > 1) {
      stats::as.formula(paste(yvar, "~ ."))
    } else {
      stats::as.formula(paste(".~", xvar))
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = negLog10p, color = factor(color_flag))) +
      ggplot2::geom_point(alpha = 0.5, size = 1.75) +
      ggplot2::geom_jitter() +
      ggplot2::labs(x = expression(log2~fold~change), y = expression(-log10~p)) +
      ggplot2::scale_color_manual(values = c(color1, color2, color3)) +
      ggplot2::facet_grid(facet_formula, space = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey80"),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "white", color = "black"),
        strip.text       = ggplot2::element_text(size = 12, face = "bold"),
        axis.title       = ggplot2::element_text(size = 12, face = "bold"),
        axis.text        = ggplot2::element_text(size = 10),
        legend.position  = "bottom"
      )
    cnt_up <- df[df$color_flag == 1, ] |>
      dplyr::group_by(.data[[xvar]], .data[[yvar]]) |> dplyr::tally(name = "n1")
    cnt_dn <- df[df$color_flag == -1, ] |>
      dplyr::group_by(.data[[xvar]], .data[[yvar]]) |> dplyr::tally(name = "n2")
    p <- p +
      ggplot2::geom_text(data = cnt_up, ggplot2::aes(label = n1), x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, color = color1) +
      ggplot2::geom_text(data = cnt_dn, ggplot2::aes(label = n2), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, color = color3)
    object@plot <- p
    invisible(object)
  }
)

#' Show a ClearScatterplot
#'
#' Prints the plot to the graphics device. If missing, builds it.
#'
#' @param object   ClearScatterplot object
#' @return         Invisibly returns the object (for chaining)
#' @export
setMethod("show", "ClearScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})


