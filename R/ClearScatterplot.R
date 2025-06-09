# Only suppress check notes on old R versions
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data", "n1", "n2"))

#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by tally

#' S4 Class: ClearScatterplot
#'
#' Encapsulates data and plotting logic for volcano plots with flexible faceting.
#'
#' @slot data   data.frame with plot data (must include `category` factor: "down"/"neutral"/"up")
#' @slot plot   ANY storing the generated ggplot object
#' @import methods
#' @exportClass ClearScatterplot
setClass(
  "ClearScatterplot",
  slots   = c(data = "data.frame", plot = "ANY"),
  validity = function(object) {
    # Required columns
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) {
      stop("Missing required columns in @data: ", paste(miss, collapse = ", "))
    }
    # We also require a semantic factor column `category`
    if (!"category" %in% names(object@data)) {
      stop("Slot `data` must contain a factor column `category` with levels 'down','neutral','up'.")
    }
    # Check exact levels
    levs <- levels(object@data$category)
    if (!identical(levs, c("down","neutral","up"))) {
      stop("`category` must have levels exactly: 'down', 'neutral', 'up'.")
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
#' Create a ClearScatterplot object from a plain data.frame of DE‐style results.
#'
#' @param data           data.frame with required columns: log2fc, negLog10p, regulation, SampleType
#' @param highLog2fc     threshold for up‐regulation (default = 0.585)
#' @param lowLog2fc      threshold for down‐regulation (default = –0.585)
#' @param negLog10pValue threshold for significance (default = 1.301)
#' @return An S4 ClearScatterplot object whose `@data` slot contains a new `category` factor
#' @export
ClearScatterplot <- function(
  data,
  highLog2fc     = 0.585,
  lowLog2fc      = -0.585,
  negLog10pValue = 1.301
) {
  stopifnot(is.data.frame(data))
  req <- c("log2fc", "negLog10p", "regulation", "SampleType")
  miss <- setdiff(req, names(data))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }
  if (!is.numeric(data$log2fc))   stop("`log2fc` must be numeric")
  if (!is.numeric(data$negLog10p)) stop("`negLog10p` must be numeric")

  na_idx <- is.na(data$log2fc) | is.na(data$negLog10p)
  if (any(na_idx)) {
    warning(sum(na_idx), " row(s) removed due to NA in `log2fc` or `negLog10p`")
    data <- data[!na_idx, , drop = FALSE]
  }

  data$color_flag <- with(
    data,
    ifelse(
      log2fc > highLog2fc & negLog10p > negLog10pValue,   1L,
      ifelse(log2fc < lowLog2fc  & negLog10p > negLog10pValue, -1L, 0L)
    )
  )
  data$category <- factor(
    data$color_flag,
    levels = c(-1L, 0L, 1L),
    labels = c("down","neutral","up"),
    ordered = TRUE
  )

  methods::new("ClearScatterplot", data = data, plot = NULL)
}

#' Internal: run DE per “cell”
#'
#' @noRd
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
    warning(sprintf("No DE results for %s / %s", cell$SampleType, cell$timePoint))
    return(NULL)
  }
  data.frame(
    log2fc     = tt$logFC,
    negLog10p  = -log10(tt$P.Value),
    regulation = ifelse(tt$logFC > 0, "up", "down"),
    SampleType = cell$SampleType,
    timePoint  = cell$timePoint,
    stringsAsFactors = FALSE,
    row.names = rownames(tt)
  )
}

#' Constructor: ClearScatterplot from MAE
#'
#' @export
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom SummarizedExperiment assay colData
#' @importFrom matrixStats rowVars
#' @importFrom stats quantile model.matrix
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

  assays <- names(MultiAssayExperiment::experiments(mae))
  if (!assayName %in% assays) {
    stop("Assay '", assayName, "' not found in this MAE.")
  }
  se   <- MultiAssayExperiment::experiments(mae)[[assayName]]
  expr <- SummarizedExperiment::assay(se)

  rv  <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
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

  shared <- intersect(colnames(expr), rownames(meta))
  if (!length(shared)) stop("No overlapping samples between assay and metadata.")
  expr <- expr[, shared, drop = FALSE]
  meta <- meta[shared, , drop = FALSE]

  keep <- rowSums(is.na(meta[needed])) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  grp_counts <- table(meta[[groupColumn]])
  if (any(grp_counts < 3)) stop("Each level of '", groupColumn, "' must have ≥ 3 samples.")
  if (ncol(expr) < 6)       stop("Too few samples (<6) remain after filtering.")

  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Constructor: ClearScatterplot_table
#'
#' @export
#' @importFrom matrixStats rowVars
#' @importFrom stats quantile
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

  rv   <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr  <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  meta2 <- meta[, needed, drop = FALSE]
  keep  <- rowSums(is.na(meta2)) == 0
  expr  <- expr[, keep, drop = FALSE]
  meta2 <- meta2[keep, , drop = FALSE]

  grp_counts <- table(meta2[[groupColumn]])
  if (any(grp_counts < 3)) stop("Each level of '", groupColumn, "' must have ≥ 3 samples.")
  if (nrow(expr) < 1)     stop("No features remain after variance filtering.")

  .ClearScatterplot_core(
    expr, meta2, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Internal core for both MAE‐ and matrix‐based constructors
#'
#' @noRd
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
    if (all(expr == floor(expr), na.rm = TRUE) && max(expr, na.rm = TRUE) > 30) {
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
      which(meta[[timepoint]]  == tp & meta[[sampleType]] == st)
    } else {
      which(meta[[sampleType]] == st)
    }
    if (length(idx) < 2 || length(unique(meta[[groupColumn]][idx])) < 2) return(NULL)
    ce <- expr[, idx, drop = FALSE]
    cm <- meta[idx, , drop = FALSE]
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) ce <- log2(ce + 1)
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) return(NULL)
    .run_DE(ce, cm, groupColumn, design, dataType, list(SampleType = st, timePoint = tp))
  }

  use_par <- (vectorized == "vectorized" ||
             (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
  lst <- if (use_par) {
    BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lapply(seq_len(nrow(cells)), run_cell)
  }

  pd <- do.call(rbind, lst)
  if (!nrow(pd)) stop("No DE results to plot (all cells returned zero rows).")
  rownames(pd) <- NULL
  ClearScatterplot(pd)
}

#' Generic: createPlot
#'
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' Method: createPlot for ClearScatterplot
#'
#' @import ggplot2
#' @exportMethod createPlot
setMethod("createPlot", "ClearScatterplot", function(object,
  color1          = "cornflowerblue",
  color2          = "grey",
  color3          = "indianred",
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

  has_tp <- "timePoint" %in% names(df) && !all(is.na(df$timePoint))
  facet <- if (has_tp && length(unique(df$timePoint)) > 1) {
    "timePoint ~ SampleType"
  } else {
    ". ~ SampleType"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = negLog10p, color = category)) +
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
         strip.text       = ggplot2::element_text(size = text_size+2, face="bold", family=text_family),
         axis.title       = ggplot2::element_text(size = text_size+2, face="bold", family=text_family),
         axis.text        = ggplot2::element_text(size = text_size, family=text_family),
         legend.position  = legend_position,
         legend.title     = ggplot2::element_text(size = text_size+2, face="bold", family=text_family),
         legend.text      = ggplot2::element_text(size = text_size, family=text_family)
       ) +
       ggplot2::scale_color_manual(
         values = c(down=color3, neutral=color2, up=color1),
         breaks = c("down","neutral","up"),
         labels = legend_labels
       )

  # add count labels
  if (has_tp && length(unique(df$timePoint))>1) {
    up <- df %>% dplyr::filter(.data$category=="up")    %>% dplyr::group_by(.data$timePoint,.data$SampleType) %>% dplyr::tally(name="n1")
    dn <- df %>% dplyr::filter(.data$category=="down")  %>% dplyr::group_by(.data$timePoint,.data$SampleType) %>% dplyr::tally(name="n2")
  } else {
    up <- df %>% dplyr::filter(.data$category=="up")    %>% dplyr::group_by(.data$SampleType) %>% dplyr::tally(name="n1")
    dn <- df %>% dplyr::filter(.data$category=="down")  %>% dplyr::group_by(.data$SampleType) %>% dplyr::tally(name="n2")
  }
  if (nrow(up)) p <- p + ggplot2::geom_text(data=up, ggplot2::aes(label=n1),
        x=Inf,y=Inf,hjust=1.1,vjust=1.1,color=color1,family=text_family,size=text_size/ggplot2::.pt)
  if (nrow(dn)) p <- p + ggplot2::geom_text(data=dn, ggplot2::aes(label=n2),
        x=-Inf,y=Inf,hjust=-0.1,vjust=1.1,color=color3,family=text_family,size=text_size/ggplot2::.pt)

  object@plot <- p
  invisible(object)
})

#' Method: show
#'
#' @exportMethod show
setMethod("show", "ClearScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})

