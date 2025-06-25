# Only suppress check notes on old R versions
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data", "n1", "n2"))

#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by tally
#' @importFrom methods setClass new setGeneric setMethod
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom BiocParallel bpparam bpworkers bplapply
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom SummarizedExperiment assay colData
#' @importFrom matrixStats rowVars
#' @importFrom stats quantile model.matrix
#' @importFrom grid unit gpar grid.points
#' @importFrom ggplot2 ggplot aes geom_jitter geom_point labs facet_grid theme_bw theme element_text element_rect scale_color_manual
#' @exportClass ThresholdedScatterplot

setClass(
  "ThresholdedScatterplot",
  slots   = c(data = "data.frame", plot = "ANY"),
  validity = function(object) {
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) {
      stop("Missing required columns in @data: ", paste(miss, collapse = ", "))
    }
    if (!"category" %in% names(object@data)) {
      stop("Slot `data` must contain factor column `category` with levels 'down','neutral','up'.")
    }
    levs <- levels(object@data$category)
    if (!identical(levs, c("down","neutral","up"))) {
      stop("`category` levels must be exactly: 'down','neutral','up'.")
    }
    TRUE
  }
)

#' @export
setGeneric("ThresholdedScatterplot", function(data, ...) standardGeneric("ThresholdedScatterplot"))

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
  if (!is.numeric(data$log2fc) || !is.numeric(data$negLog10p)) {
    stop("`log2fc` and `negLog10p` must be numeric")
  }

  na_idx <- is.na(data$log2fc) | is.na(data$negLog10p)
  if (any(na_idx)) {
    warning(sum(na_idx), " row(s) removed due to NA")
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

  methods::new("ThresholdedScatterplot", data = data, plot = NULL)
}

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
  if (!assayName %in% names(assays)) {
    stop("Assay '", assayName, "' not found in this MAE.")
  }
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
  if (any(grp_counts < 3)) stop(
    sprintf("Each level of '%s' must have ≥ 3 samples.", groupColumn)
  )
  if (ncol(expr) < 6) stop("Too few samples (<6) remain after filtering.")

  .ThresholdedScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

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
  if (any(grp_counts < 3)) {
    stop(sprintf("Each level of '%s' must have ≥ 3 samples.", groupColumn))
  }
  if (nrow(expr) < 1) {
    stop("No features remain after variance filtering.")
  }

  .ThresholdedScatterplot_core(
    expr, meta2, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#’ @title Faceted Volcano from Multiple Inputs (max flexibility)
#’ @description
#’ Compute and plot volcano scatterplots for multiple conditions in one faceted ggplot.  
#’ Supports three input types:
#’   - “expr”: named list of expression matrices + matching metadata frames  
#’   - “mae”:  named list of MultiAssayExperiment objects  
#’   - “de”:   named list of precomputed DE tables  
#’
#’ You may pass extra parameters to the DE‐calculation step via `compute_args`
#’ and to the plotting step via `plot_args`.
#’
#’ @param data_list    Named list of inputs (matrices, MAEs, or DE tables).
#’ @param meta_list    Named list of metadata data.frames (required for expr mode).
#’ @param input_type   One of `"expr"`, `"mae"`, or `"de"`.
#’ @param groupColumn  Grouping column (expr/mae only).
#’ @param sampleType   SampleType column (expr/mae only).
#’ @param timepoint    Optional timepoint column (expr/mae only).
#’ @param dataType     `"auto"`, `"continuous"` or `"count"` (expr/mae only).
#’ @param var_quantile Variance quantile in [0,1] (expr/mae only).
#’ @param compute_args Named list of extra args passed to the DE function.
#’ @param plot_args    Named list of extra args passed to `createPlot()`.
#’                    (*Supported:* color1, color2, color3, point_size, point_alpha, text_size, etc.)
#’ @param facet_by     Facet formula, e.g. `". ~ panel"` or `"panel ~ SampleType"`.
#’ @return A `ThresholdedScatterplot` S4 object whose `@plot` slot is the faceted volcano.
#’ @importFrom ggplot2 facet_grid geom_text
#’ @importFrom dplyr filter group_by tally
#’ @importFrom rlang %||%
#’ @importFrom MultiAssayExperiment experiments
#’ @export
ThresholdedScatterplot_list <- function(
  data_list,
  meta_list     = NULL,
  input_type    = c("expr","mae","de"),
  groupColumn   = "Group",
  sampleType    = "SampleType",
  timepoint     = NULL,
  dataType      = c("auto","continuous","count"),
  var_quantile  = 0.75,
  var_quantile     = 0.75,
  highLog2fc       = 0.585,
  lowLog2fc        = -0.585,
  negLog10pValue   = 1.301,
  compute_args  = list(),
  plot_args     = list(),
  facet_by      = ". ~ panel"
) {
  ## 1) validate
  input_type <- match.arg(input_type)
  dataType   <- match.arg(dataType)
  stopifnot(is.list(data_list), length(data_list)>0, !is.null(names(data_list)))
  if (length(compute_args) && is.null(names(compute_args)))
    stop("`compute_args` must be a named list.")
  if (length(plot_args) && is.null(names(plot_args)))
    stop("`plot_args` must be a named list.")
  panels <- names(data_list)
  if (input_type=="expr") {
    stopifnot(is.list(meta_list), identical(sort(names(meta_list)), sort(panels)))
    lapply(meta_list, function(df) stopifnot(is.data.frame(df)))
  }

  ## 2) helper: get one DE table
  get_de <- function(x, meta) {
    if (input_type == "expr") {
      stopifnot(is.matrix(x), is.data.frame(meta))
      args <- c(list(
        expr        = x,
        meta        = meta,
        groupColumn = groupColumn,
        sampleType  = sampleType,
        timepoint   = timepoint,
        dataType    = dataType,
        var_quantile= var_quantile
      ), compute_args)
      tbl <- do.call(ThresholdedScatterplot_table, args)@data

    } else if (input_type == "mae") {
      stopifnot(inherits(x, "MultiAssayExperiment"))
      assays    <- MultiAssayExperiment::experiments(x)
      assayName <- names(assays)[1]
      args <- c(list(
        mae         = x,
        assayName   = assayName,
        groupColumn = groupColumn,
        sampleType  = sampleType,
        timepoint   = timepoint,
        dataType    = dataType,
        var_quantile= var_quantile
      ), compute_args)
      tbl <- do.call(ThresholdedScatterplot_MAE, args)@data

    } else {  # "de"
      stopifnot(is.data.frame(x))
      tbl <- x
      req <- c("log2fc","negLog10p","regulation","SampleType")
      miss <- setdiff(req, colnames(tbl))
      if (length(miss)) stop("Missing DE columns: ", paste(miss, collapse=", "))
    }
    if (nrow(tbl)==0) warning("Panel yields zero DE rows; it will be dropped.")
    tbl
  }

  ## 3) build & drop empty panels (use Map so names stay aligned)
  de_list <- Map(function(name, x) {
    md <- if (input_type=="expr") meta_list[[name]] else NULL
    df <- get_de(x, md)
    if (nrow(df)==0) return(NULL)
    df$panel <- name
    df
  }, panels, data_list)
  de_list <- Filter(Negate(is.null), de_list)
  if (length(de_list)==0) stop("All panels dropped; no DE results to plot.")
  all_de <- do.call(rbind, de_list)
  all_de$panel <- factor(all_de$panel, levels = names(de_list))

  ## 4) initial S4 object + ggplot
  cs      <- ThresholdedScatterplot(
                                  all_de,
                                  highLog2fc = highLog2fc,
                                  lowLog2fc = lowLog2fc,
                                  negLog10pValue = negLog10pValue
                                  )
           
  cs_full <- do.call(createPlot, c(list(cs), plot_args))
  p       <- cs_full@plot

  ## 5) strip out the built-in n1/n2 labels
  p$layers <- Filter(function(ly) !any(c("n1","n2") %in% names(ly$data)), p$layers)

  ## 6) recompute per-panel up/down counts
  counts_up <- all_de %>%
    dplyr::filter(regulation=="up") %>%
    dplyr::group_by(panel, SampleType, timePoint) %>%
    dplyr::tally(name="n_up")
  counts_dn <- all_de %>%
    dplyr::filter(regulation=="down") %>%
    dplyr::group_by(panel, SampleType, timePoint) %>%
    dplyr::tally(name="n_dn")

  ## 7) extract user-colors (or defaults)
  col_up   <- plot_args$color1    %||% "cornflowerblue"
  col_dn   <- plot_args$color3    %||% "indianred"
  txt_sz   <- plot_args$text_size %||% 10

  ## 8) add correct per-panel count labels
  p <- p +
    ggplot2::geom_text(
      data = counts_up,
      ggplot2::aes(x=Inf,  y=Inf, label=n_up),
      hjust=1.1, vjust=1.1, color=col_up, size=txt_sz/ggplot2::.pt
    ) +
    ggplot2::geom_text(
      data = counts_dn,
      ggplot2::aes(x=-Inf, y=Inf, label=n_dn),
      hjust=-0.1, vjust=1.1, color=col_dn, size=txt_sz/ggplot2::.pt
    )

  ## 9) override facet with your `facet_by`
  facet_formula <- tryCatch(
    stats::as.formula(facet_by),
    error = function(e) stop("`facet_by` is not a valid formula: ", facet_by)
  )
  cs@plot <- p + ggplot2::facet_grid(facet_formula, space="free", scales="free_x")

  invisible(cs)
}


           
#' @noRd
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
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) {
      ce <- log2(ce + 1)
    }
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) return(NULL)
    .run_DE(ce, cm, groupColumn, design, dataType, list(
      SampleType = st, timePoint = tp
    ))
  }

  use_par <- (vectorized == "vectorized" ||
             (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
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

#' @exportMethod createPlot
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' @describeIn createPlot Render the volcano ggplot
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
    color = color1, family = text_family, size = text_size/ggplot2::.pt
  )
  if (nrow(dn)) p <- p + ggplot2::geom_text(
    data = dn, ggplot2::aes(label = n2),
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
    color = color3, family = text_family, size = text_size/ggplot2::.pt
  )

  object@plot <- p
  invisible(object)
})

#' @exportMethod show
setMethod("show", "ThresholdedScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})
