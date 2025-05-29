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
#'
#' @param data data.frame
#' @param ... additional arguments
#' @export
setGeneric("ClearScatterplot", function(data, ...) standardGeneric("ClearScatterplot"))

#' Constructor: ClearScatterplot
#'
#' Create a ClearScatterplot S4 object from a precomputed data frame. 
#' Automatically flags up/down-regulated genes.
#'
#' @param data           data.frame with required columns
#' @param highLog2fc     numeric threshold for up-regulation (default 0.585)
#' @param lowLog2fc      numeric threshold for down-regulation (default -0.585)
#' @param negLog10pValue numeric threshold for p-value significance (default 1.301)
#' @return ClearScatterplot instance
#' @examples
#' df <- data.frame(
#'   log2fc = rnorm(100),
#'   negLog10p = runif(100),
#'   regulation = sample(c("up","down"),100,TRUE),
#'   SampleType = sample(c("A","B"),100,TRUE),
#'   timePoint = sample(c("T1","T2"),100,TRUE)
#' )
#' cs <- ClearScatterplot(df)
#' createPlot(cs)
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
#'
#' @keywords internal
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
#' This function builds a ClearScatterplot object from a MultiAssayExperiment, 
#' running limma DE per SampleType x TimePoint cell (if present).
#'
#' @param mae         MultiAssayExperiment
#' @param assayName   character: name of assay
#' @param groupColumn character: column in colData defining groups
#' @param sampleType  character: column for faceting (columns)
#' @param timepoint   character or NULL: column for faceting (rows)
#' @param dataType    one of "auto", "continuous", "count"
#' @param vectorized  one of "auto", "perCell", "vectorized"
#' @param BPPARAM     BiocParallelParam for parallelization
#' @return ClearScatterplot
#'
#' @examples
#' library(MultiAssayExperiment)
#' data("miniACC", package="MultiAssayExperiment")
#'
#' # Remove NAs from key metadata columns before analysis
#' md <- as.data.frame(colData(miniACC)[,c("SCNA.cluster","pathologic_stage","C1A.C1B")])
#' keep <- rownames(md)[!is.na(md$SCNA.cluster) & !is.na(md$pathologic_stage) & !is.na(md$C1A.C1B)]
#' miniACC_clean <- miniACC[, keep]
#'
#' cs <- ClearScatterplot_MAE(
#'   mae         = miniACC_clean,
#'   assayName   = "RNASeq2GeneNorm",
#'   sampleType  = "pathologic_stage",
#'   timepoint   = "SCNA.cluster",
#'   groupColumn = "C1A.C1B",
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
    expr       = expr,
    meta       = meta,
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
#' This function builds a ClearScatterplot object from a numeric matrix 
#' (e.g. normalized gene expression or counts) and metadata table, 
#' running limma DE per SampleType x TimePoint cell (if present).
#'
#' @param expr        matrix (rows: features, cols: samples)
#' @param meta        data.frame (metadata, rownames should match colnames(expr))
#' @param groupColumn character: group variable
#' @param sampleType  character: facet columns
#' @param timepoint   character or NULL: facet rows
#' @param dataType    one of "auto", "continuous", "count"
#' @param vectorized  one of "auto", "perCell", "vectorized"
#' @param BPPARAM     BiocParallelParam
#' @return ClearScatterplot
#'
#' @examples
#' # Example: direct from table
#' expr <- matrix(rpois(1000, lambda=10), nrow=100)
#' meta <- data.frame(
#'   SampleType = rep(c("Tissue","Blood"), each=5),
#'   timePoint = rep(c("T1","T2"), 5),
#'   Group = rep(c("Case","Control"), 5),
#'   row.names = paste0("S", 1:10)
#' )
#' cs <- ClearScatterplot_table(
#'   expr = expr[,1:10],
#'   meta = meta[1:10,],
#'   groupColumn = "Group",
#'   sampleType = "SampleType",
#'   timepoint = "timePoint"
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
    expr       = expr,
    meta       = meta,
    groupColumn = groupColumn,
    sampleType  = sampleType,
    timepoint   = timepoint,
    dataType    = dataType,
    vectorized  = vectorized,
    BPPARAM     = BPPARAM
  )
}

#' Internal core: shared between MAE and table constructors
#'
#' @keywords internal
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

#' S4 Generic: createPlot
#'
#' @param object ClearScatterplot object
#' @param ... extra arguments (see method)
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' S4 Method: createPlot for ClearScatterplot
#'
#' Render the volcano scatterplot, faceted by SampleType and timePoint if present, with up/down gene count overlays.
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





