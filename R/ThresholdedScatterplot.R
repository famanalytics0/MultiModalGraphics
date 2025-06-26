# Only suppress check notes on old R versions
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data", "n1", "n2"))

# Null-coalesce operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by tally filter
#' @importFrom methods setClass new setGeneric setMethod
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom BiocParallel bpparam bpworkers bplapply
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom SummarizedExperiment assay colData
#' @importFrom matrixStats rowVars
#' @importFrom stats quantile model.matrix
#' @importFrom ggplot2 ggplot aes geom_jitter geom_point geom_text labs facet_grid theme_bw theme element_text element_rect scale_color_manual

#' @title S4 Class for Thresholded Volcano Scatterplots
#' @description S4 class and methods for thresholded volcano scatterplots of DE results.
#' @slot data Data frame with required columns: log2fc, negLog10p, regulation, SampleType, category.
#' @slot plot ggplot2 object, NULL if not yet created.
#' @exportClass ThresholdedScatterplot
setClass(
  "ThresholdedScatterplot",
  slots = c(data = "data.frame", plot = "ANY"),
  validity = function(object) {
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) stop("Missing required columns in @data: ", paste(miss, collapse = ", "))
    if (!"category" %in% names(object@data)) stop("Slot `data` must contain factor column `category` with levels 'down','neutral','up'.")
    if (any(is.na(object@data$category)))
      stop("Column `category` contains NA values. Check input or thresholds.")
    levs <- levels(object@data$category)
    if (!identical(levs, c("down", "neutral", "up"))) stop("`category` levels must be exactly: 'down','neutral','up'.")
    if (!is.character(object@data$regulation) && !is.factor(object@data$regulation))
      stop("Column `regulation` must be character or factor.")
    TRUE
  }
)

#' @export
setGeneric("ThresholdedScatterplot", function(data, ...) standardGeneric("ThresholdedScatterplot"))

#' Accessor for data slot
#' @export
setGeneric("getData", function(object) standardGeneric("getData"))
#' @describeIn ThresholdedScatterplot Get data slot
setMethod("getData", "ThresholdedScatterplot", function(object) object@data)

#' Accessor for plot slot
#' @export
setGeneric("getPlot", function(object) standardGeneric("getPlot"))
#' @describeIn ThresholdedScatterplot Get plot slot
setMethod("getPlot", "ThresholdedScatterplot", function(object) object@plot)

#' Setter for data slot (useful for programmatic update)
#' @export
setGeneric("setData", function(object, value) standardGeneric("setData"))
#' @describeIn ThresholdedScatterplot Set data slot
setMethod("setData", "ThresholdedScatterplot", function(object, value) {
  object@data <- value
  validObject(object)
  object
})

#' Setter for plot slot (useful for programmatic update)
#' @export
setGeneric("setPlot", function(object, value) standardGeneric("setPlot"))
#' @describeIn ThresholdedScatterplot Set plot slot
setMethod("setPlot", "ThresholdedScatterplot", function(object, value) {
  object@plot <- value
  object
})

#' Create a ThresholdedScatterplot object from DE data
#'
#' @param data Data.frame with required columns: log2fc, negLog10p, regulation, SampleType
#' @param highLog2fc,lowLog2fc,negLog10pValue Numeric cutoffs for thresholding
#' @return S4 object of class ThresholdedScatterplot
#' @export
ThresholdedScatterplot <- function(
  data,
  highLog2fc     = 0.585,
  lowLog2fc      = -0.585,
  negLog10pValue = 1.301
) {
  stopifnot(is.data.frame(data))
  req <- c("log2fc", "negLog10p", "regulation", "SampleType")
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  if (!is.numeric(data$log2fc) || !is.numeric(data$negLog10p)) stop("`log2fc` and `negLog10p` must be numeric")
  if (!is.character(data$regulation) && !is.factor(data$regulation))
    stop("`regulation` column must be character or factor")
  na_idx <- is.na(data$log2fc) | is.na(data$negLog10p)
  if (any(na_idx)) {
    warning(sum(na_idx), " row(s) removed due to NA")
    data <- data[!na_idx, , drop = FALSE]
  }
  data$color_flag <- with(
    data,
    ifelse(
      log2fc > highLog2fc & negLog10p > negLog10pValue,   1L,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1L, 0L)
    )
  )
  data$category <- factor(
    data$color_flag,
    levels = c(-1L, 0L, 1L),
    labels = c("down", "neutral", "up"),
    ordered = TRUE
  )
  methods::new("ThresholdedScatterplot", data = data, plot = NULL)
}

#' @exportMethod createPlot
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' Create the volcano plot for a ThresholdedScatterplot object
#'
#' @param object ThresholdedScatterplot S4 object
#' @param color1, color2, color3 Colors for 'up', 'neutral', 'down'
#' @param xlab, ylab Axis labels
#' @param custom_theme Additional ggplot2 theme() or element
#' @param point_alpha, point_size Transparency/size for points
#' @param jitter Enable/disable jitter layer
#' @param geom_point Enable/disable plain point layer
#' @param legend_position, legend_title, legend_labels, text_family, text_size Legend and text styling
#' @param ... Passed to ggplot2::labs()
#' @return The object (invisible)
#' @exportMethod createPlot
setMethod("createPlot", "ThresholdedScatterplot", function(object,
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
  jitter          = TRUE,
  geom_point      = TRUE,
  ...
) {
  df <- object@data
  # Defensive: create timePoint column if missing
  if (!"timePoint" %in% colnames(df)) df$timePoint <- NA_character_
  has_tp <- !all(is.na(df$timePoint))
  facet <- if (has_tp && length(unique(df$timePoint)) > 1) {
    "timePoint ~ SampleType"
  } else {
    ". ~ SampleType"
  }
  if (!is.null(legend_labels) && length(legend_labels) != 3) stop("legend_labels must be NULL or length 3.")
  .pt <- tryCatch(get(".pt", envir = asNamespace("ggplot2")), error = function(e) 2.845276)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = negLog10p, color = category))
  if (isTRUE(jitter))    p <- p + ggplot2::geom_jitter(alpha = point_alpha, size = point_size)
  if (isTRUE(geom_point)) p <- p + ggplot2::geom_point(alpha = point_alpha, size = point_size)
  p <- p +
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
      values = c(down = color3, neutral = color2, up = color1),
      breaks = c("down", "neutral", "up"),
      labels = legend_labels
    )
  if (!is.null(custom_theme)) p <- p + custom_theme
  up <- df %>%
    dplyr::filter(.data$category == "up") %>%
    dplyr::group_by(.data$SampleType, .data$timePoint) %>%
    dplyr::tally(name = "n1")
  dn <- df %>%
    dplyr::filter(.data$category == "down") %>%
    dplyr::group_by(.data$SampleType, .data$timePoint) %>%
    dplyr::tally(name = "n2")
  if (nrow(up)) p <- p + ggplot2::geom_text(
    data = up, ggplot2::aes(label = n1),
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
    color = color1, family = text_family, size = text_size/.pt
  )
  if (nrow(dn)) p <- p + ggplot2::geom_text(
    data = dn, ggplot2::aes(label = n2),
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
    color = color3, family = text_family, size = text_size/.pt
  )
  object@plot <- p
  invisible(object)
})

#' @exportMethod show
#' @describeIn ThresholdedScatterplot Show method for S4 object
setMethod("show", "ThresholdedScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})

#' @title Internal: Run DE per cell
#' @param expr matrix
#' @param meta data.frame
#' @param groupColumn column for group
#' @param design design matrix
#' @param dataType 'continuous' or 'count'
#' @param cell list with SampleType, timePoint
#' @return data.frame of DE results
#' @keywords internal
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

#' Flexible table input for ThresholdedScatterplot
#'
#' @param expr Numeric matrix (rows = features, columns = samples)
#' @param meta Data.frame with sample annotations
#' @param groupColumn Grouping column name
#' @param sampleType Sample type column name
#' @param timepoint Optional timepoint column name
#' @param dataType Data type for DE ("auto","continuous","count")
#' @param vectorized Use vectorized (parallel) DE
#' @param parallel Whether to use parallel backend
#' @param BPPARAM BiocParallel param
#' @param var_quantile Row variance quantile for feature filtering
#' @return S4 object of class ThresholdedScatterplot
#' @export
ThresholdedScatterplot_table <- function(
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
  if (any(grp_counts < 3)) stop(sprintf("Each level of '%s' must have ≥ 3 samples.", groupColumn))
  if (nrow(expr) < 1) stop("No features remain after variance filtering.")
  .ThresholdedScatterplot_core(
    expr, meta2, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Flexible MultiAssayExperiment input for ThresholdedScatterplot
#'
#' @param mae MultiAssayExperiment
#' @param assayName Name of assay
#' @param groupColumn Grouping column
#' @param sampleType Sample type column
#' @param timepoint Timepoint column
#' @param dataType Data type
#' @param vectorized Vectorized DE
#' @param parallel Use parallel
#' @param BPPARAM BiocParallel param
#' @param var_quantile Row variance quantile for filtering
#' @return S4 object of class ThresholdedScatterplot
#' @export
ThresholdedScatterplot_MAE <- function(
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
  assays <- MultiAssayExperiment::experiments(mae)
  if (!assayName %in% names(assays)) stop("Assay '", assayName, "' not found in this MAE.")
  se   <- assays[[assayName]]
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
  if (any(grp_counts < 3)) stop(sprintf("Each level of '%s' must have ≥ 3 samples.", groupColumn))
  if (ncol(expr) < 6) stop("Too few samples (<6) remain after filtering.")
  .ThresholdedScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' @title Faceted Volcano from Multiple Inputs (max flexibility)
#' @param data_list List of matrix/data.frame/MultiAssayExperiment objects
#' @param meta_list List of meta data.frames (if needed)
#' @param input_type Type of input: "expr", "mae", "de"
#' @param groupColumn ...
#' @param sampleType ...
#' @param timepoint ...
#' @param dataType ...
#' @param var_quantile ...
#' @param highLog2fc ...
#' @param lowLog2fc ...
#' @param negLog10pValue ...
#' @param compute_args ...
#' @param plot_args ...
#' @param facet_by ...
#' @return S4 object of class ThresholdedScatterplot
#' @export
ThresholdedScatterplot_list <- function(
  data_list,
  meta_list       = NULL,
  input_type      = c("expr", "mae", "de"),
  groupColumn     = "Group",
  sampleType      = "SampleType",
  timepoint       = NULL,
  dataType        = c("auto", "continuous", "count"),
  var_quantile    = 0.75,
  highLog2fc      = 0.585,
  lowLog2fc       = -0.585,
  negLog10pValue  = 1.301,
  compute_args    = list(),
  plot_args       = list(),
  facet_by        = ". ~ panel"
) {
  input_type <- match.arg(input_type)
  dataType   <- match.arg(dataType)
  stopifnot(is.list(data_list), length(data_list) > 0, !is.null(names(data_list)))
  stopifnot(is.list(compute_args), is.list(plot_args))
  if (length(compute_args) && is.null(names(compute_args))) stop("`compute_args` must be a named list.")
  if (length(plot_args) && is.null(names(plot_args))) stop("`plot_args` must be a named list.")
  panels <- names(data_list)
  if (input_type == "expr") {
    stopifnot(is.list(meta_list), identical(sort(names(meta_list)), sort(panels)))
    lapply(meta_list, function(df) stopifnot(is.data.frame(df)))
  }
  get_de <- function(x, meta) {
    if (input_type == "expr") {
      stopifnot(is.matrix(x), is.data.frame(meta))
      args <- c(list(expr = x, meta = meta, groupColumn = groupColumn, sampleType  = sampleType,
                     timepoint   = timepoint, dataType = dataType, var_quantile= var_quantile), compute_args)
      tbl <- do.call(ThresholdedScatterplot_table, args)@data
    } else if (input_type == "mae") {
      stopifnot(inherits(x, "MultiAssayExperiment"))
      assays    <- MultiAssayExperiment::experiments(x)
      assayName <- names(assays)[1]
      args <- c(list(mae = x, assayName = assayName, groupColumn = groupColumn, sampleType  = sampleType,
                     timepoint   = timepoint, dataType = dataType, var_quantile= var_quantile), compute_args)
      tbl <- do.call(ThresholdedScatterplot_MAE, args)@data
    } else {
      stopifnot(is.data.frame(x))
      tbl <- x
      req <- c("log2fc","negLog10p","regulation","SampleType","category")
      miss <- setdiff(req, colnames(tbl))
      if (length(miss)) stop("Missing DE columns: ", paste(miss, collapse=", "))
    }
    if (nrow(tbl) == 0) warning("Panel yields zero DE rows; it will be dropped.")
    tbl
  }
  de_list <- Map(function(name, x) {
    md <- if (input_type == "expr") meta_list[[name]] else NULL
    df <- get_de(x, md)
    if (nrow(df) == 0) return(NULL)
    df$panel <- name
    df
  }, panels, data_list)
  de_list <- Filter(Negate(is.null), de_list)
  if (length(de_list) == 0) stop("All panels dropped; no DE results to plot.")
  all_de <- do.call(rbind, de_list)
  all_de$panel <- factor(all_de$panel, levels = names(de_list))
  if (!"category" %in% colnames(all_de))
    stop("Input DE tables must include a 'category' column. Did you use the correct DE wrapper?")
  if (!is.factor(all_de$category))
    all_de$category <- factor(all_de$category, levels = c("down", "neutral", "up"), ordered = TRUE)
  tp_col <- grep("^time[Pp]oint$", colnames(all_de), value = TRUE)
  if (is.null(timepoint) && length(tp_col) == 1) {
    all_de[[tp_col]] <- NULL
  }
  if (!is.null(timepoint) && length(tp_col) != 1)
    stop("Could not find a unique 'timepoint' column for overlays.")
  cs <- ThresholdedScatterplot(
    all_de,
    highLog2fc     = highLog2fc,
    lowLog2fc      = lowLog2fc,
    negLog10pValue = negLog10pValue
  )
  cs_full <- do.call(createPlot, c(list(cs), plot_args))
  p       <- cs_full@plot
  p$layers <- Filter(function(ly) {
    dat <- tryCatch(ly$data, error = function(e) NULL)
    !(is.data.frame(dat) && any(c("n1", "n2") %in% names(dat)))
  }, p$layers)
  if (is.null(timepoint)) {
    counts_up <- all_de %>%
      dplyr::filter(category == "up") %>%
      dplyr::group_by(panel, SampleType) %>%
      dplyr::tally(name = "n_up")
    counts_dn <- all_de %>%
      dplyr::filter(category == "down") %>%
      dplyr::group_by(panel, SampleType) %>%
      dplyr::tally(name = "n_dn")
  } else {
    counts_up <- all_de %>%
      dplyr::filter(category == "up") %>%
      dplyr::group_by(panel, SampleType, .data[[tp_col]]) %>%
      dplyr::tally(name = "n_up")
    counts_dn <- all_de %>%
      dplyr::filter(category == "down") %>%
      dplyr::group_by(panel, SampleType, .data[[tp_col]]) %>%
      dplyr::tally(name = "n_dn")
  }
  col_up <- plot_args$color1    %||% "cornflowerblue"
  col_dn <- plot_args$color3    %||% "indianred"
  txt_sz <- plot_args$text_size %||% 10L
  .pt <- tryCatch(get(".pt", envir = asNamespace("ggplot2")), error = function(e) 2.845276)
  p <- p +
    ggplot2::geom_text(
      data = counts_up,
      ggplot2::aes(x = Inf, y = Inf, label = n_up),
      hjust = 1.1, vjust = 1.1, color = col_up, size = txt_sz / .pt
    ) +
    ggplot2::geom_text(
      data = counts_dn,
      ggplot2::aes(x = -Inf, y = Inf, label = n_dn),
      hjust = -0.1, vjust = 1.1, color = col_dn, size = txt_sz / .pt
    )
  facet_formula <- tryCatch(
    stats::as.formula(facet_by),
    error = function(e) stop("`facet_by` is not a valid formula: ", facet_by)
  )
  cs@plot <- p + ggplot2::facet_grid(facet_formula, space = "free", scales = "free_x")
  invisible(cs)
}

#' Internal: run DE for each cell of stratification
#' @param expr Expression matrix
#' @param meta Meta data
#' @param groupColumn Group column
#' @param sampleType Sample type
#' @param timepoint Time point
#' @param dataType Data type
#' @param vectorized Vectorized flag
#' @param BPPARAM BiocParallel param
#' @return S4 object of class ThresholdedScatterplot
#' @keywords internal
.ThresholdedScatterplot_core <- function(
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
      which(meta[[timepoint]] == tp & meta[[sampleType]] == st)
    } else {
      which(meta[[sampleType]] == st)
    }
    if (length(idx) < 2 || length(unique(meta[[groupColumn]][idx])) < 2) return(NULL)
    ce <- expr[, idx, drop = FALSE]
    cm <- meta[idx, , drop = FALSE]
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) ce <- log2(ce + 1)
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) return(NULL)
    .run_DE(ce, cm, groupColumn, design, dataType, list(
      SampleType = st, timePoint = tp
    ))
  }
  use_par <- (vectorized == "vectorized" || (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
  lst <- if (use_par) {
    BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lapply(seq_len(nrow(cells)), run_cell)
  }
  pd <- do.call(rbind, lst)
  if (is.null(pd) || nrow(pd) == 0) stop("No DE results to plot (all cells returned zero rows).")
  rownames(pd) <- NULL
  ThresholdedScatterplot(pd)
}

