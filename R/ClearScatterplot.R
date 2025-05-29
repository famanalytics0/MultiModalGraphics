# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2"))

#' S4 Class: ClearScatterplot
#'
#' Encapsulates data and plotting logic for volcano plots with flexible faceting.
#'
#' @slot data data.frame with columns: log2fc, negLog10p, regulation, SampleType, (optional) timePoint, color_flag
#' @slot plot ANY storing the generated ggplot object
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

#' S4 Generic: ClearScatterplot
#' @param data data.frame
#' @param ... additional arguments
#' @export
setGeneric("ClearScatterplot", function(data, ...) standardGeneric("ClearScatterplot"))

#' Constructor: ClearScatterplot
#'
#' @param data           data.frame with required columns
#' @param highLog2fc     numeric threshold for up-regulation (default 0.585)
#' @param lowLog2fc      numeric threshold for down-regulation (default -0.585)
#' @param negLog10pValue numeric threshold for p-value significance (default 1.301)
#' @return ClearScatterplot instance
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
  methods::new("ClearScatterplot", data = data, plot = NULL)
}

#' Internal utility: Differential Expression for a cell (SampleType x TimePoint)
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

#' Constructor: ClearScatterplot from MultiAssayExperiment
#'
#' This function does all NA removal, matching, and sample sanity checks internally,
#' exactly as specified by the workflow you requested (feature filter, NA filter, etc).
#'
#' @param mae         MultiAssayExperiment
#' @param assayName   character
#' @param groupColumn character
#' @param sampleType  character
#' @param timepoint   character or NULL
#' @param dataType    "auto", "continuous", "count"
#' @param vectorized  "auto", "perCell", "vectorized"
#' @param varFilterQuantile numeric. Quantile for row variance filtering (default 0.75)
#' @param BPPARAM     BiocParallelParam
#' @return ClearScatterplot
#' @export
ClearScatterplot_MAE <- function(
  mae,
  assayName,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto", "continuous", "count"),
  vectorized  = c("auto", "perCell", "vectorized"),
  varFilterQuantile = 0.75,
  BPPARAM     = BiocParallel::bpparam()
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)

  # 1. Pull expression matrix and feature filter
  expr_full <- SummarizedExperiment::assay(
    MultiAssayExperiment::experiments(mae)[[assayName]]
  )
  rowvars <- matrixStats::rowVars(expr_full)
  keep_feat <- which(rowvars >= quantile(rowvars, varFilterQuantile, na.rm = TRUE))
  expr_filt <- expr_full[keep_feat, , drop = FALSE]

  # 2. Subset MAE to filtered features
  mae_filt <- subsetByRow(mae, S4Vectors::IntegerList(setNames(list(keep_feat), assayName)))

  # 3. Pull meta for relevant columns only
  meta_cols <- unique(c(groupColumn, sampleType, timepoint))
  meta <- as.data.frame(
    SummarizedExperiment::colData(mae_filt)[, meta_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  # 4. Vectorized: retain only samples with no NA in needed columns
  keep_samps <- rownames(meta)[rowSums(is.na(meta)) == 0]
  if(length(keep_samps) < 6) stop("Too few samples after NA removal.")
  mae_clean  <- mae_filt[, keep_samps]
  expr_clean <- SummarizedExperiment::assay(
    MultiAssayExperiment::experiments(mae_clean)[[assayName]]
  )
  meta_clean <- as.data.frame(
    SummarizedExperiment::colData(mae_clean)[, meta_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  # 5. Group size check
  n_group <- table(meta_clean[[groupColumn]])
  if(any(n_group < 3)) stop("Fewer than 3 samples per group after NA removal.")

  .ClearScatterplot_core(
    expr       = expr_clean,
    meta       = meta_clean,
    groupColumn = groupColumn,
    sampleType  = sampleType,
    timepoint   = timepoint,
    dataType    = dataType,
    vectorized  = vectorized,
    BPPARAM     = BPPARAM
  )
}

#' Constructor: ClearScatterplot from Expression Matrix + Metadata Table
#'
#' Handles NA/matching inside.
#'
#' @param expr        matrix
#' @param meta        data.frame
#' @param groupColumn character
#' @param sampleType  character
#' @param timepoint   character or NULL
#' @param dataType    "auto", "continuous", "count"
#' @param vectorized  "auto", "perCell", "vectorized"
#' @param BPPARAM     BiocParallelParam
#' @return ClearScatterplot
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

  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  meta_cols <- needed

  # Only samples in both meta and expr
  shared_samples <- intersect(rownames(meta), colnames(expr))
  meta <- meta[shared_samples, , drop = FALSE]
  expr <- expr[, shared_samples, drop = FALSE]

  # Apply the same fully vectorized NA-removal
  keep_samps <- rownames(meta)[rowSums(is.na(meta[, meta_cols, drop = FALSE])) == 0]
  if(length(keep_samps) < 6) stop("Too few samples after NA removal.")
  expr_clean <- expr[, keep_samps, drop = FALSE]
  meta_clean <- meta[keep_samps, , drop = FALSE]

  # Group size check
  n_group <- table(meta_clean[[groupColumn]])
  if(any(n_group < 3)) stop("Fewer than 3 samples per group after NA removal.")

  .ClearScatterplot_core(
    expr       = expr_clean,
    meta       = meta_clean,
    groupColumn = groupColumn,
    sampleType  = sampleType,
    timepoint   = timepoint,
    dataType    = dataType,
    vectorized  = vectorized,
    BPPARAM     = BPPARAM
  )
}

#' Internal core: shared between MAE and table constructors
.ClearScatterplot_core <- function(
  expr, meta, groupColumn, sampleType, timepoint,
  dataType, vectorized, BPPARAM
) {
  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  keep <- rownames(meta)[rowSums(is.na(meta[, needed, drop = FALSE])) == 0]
  if (length(keep) < 6) stop("Too few samples after NA/matching removal.")
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

#' S4 Generic: createPlot
#'
#' @param object ClearScatterplot object
#' @param ... extra arguments (see method)
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' S4 Method: createPlot for ClearScatterplot
#'
#' @param object   ClearScatterplot
#' @param color1   up-regulation color
#' @param color2   neutral color
#' @param color3   down-regulation color
#' @return invisibly returns object with plot slot filled
#' @exportMethod createPlot
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

#' S4 Method: show for ClearScatterplot
#'
#' Prints the volcano plot.
#' @param object ClearScatterplot
#' @exportMethod show
setMethod("show", "ClearScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})





