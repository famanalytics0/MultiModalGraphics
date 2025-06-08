# Only suppress check notes on old R versions
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data", "n1", "n2"))

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
  slots = c(data = "data.frame", plot = "ANY"),
  validity = function(object) {
    # Required columns
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) {
      stop("Missing required columns in @data: ", paste(miss, collapse = ", "))
    }
    # We also require a semantic factor column `category`
    if (!"category" %in% names(object@data)) {
      stop("Slot `data` must contain a factor column `category` with levels 'down', 'neutral', 'up'.")
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
  # 1) Check that the required columns are present before doing any type checks
  stopifnot(is.data.frame(data))
  req <- c("log2fc", "negLog10p", "regulation", "SampleType")
  miss <- setdiff(req, names(data))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }
  # 2) Now that columns exist, confirm types
  if (!is.numeric(data$log2fc))    stop("`log2fc` must be numeric")
  if (!is.numeric(data$negLog10p))  stop("`negLog10p` must be numeric")
  # (We do not strictly type‐check `regulation` or `SampleType` here;
  #  they can be character or factor, as long as they exist.)

  # 3) Remove rows with NA in key numeric columns
  na_idx <- is.na(data$log2fc) | is.na(data$negLog10p)
  if (any(na_idx)) {
    warning(sum(na_idx), " row(s) removed due to NA in `log2fc` or `negLog10p`")
    data <- data[!na_idx, , drop = FALSE]
  }

  # 4) Compute a numeric “color_flag” for threshold‐based categorization
  data$color_flag <- with(
    data,
    ifelse(
      log2fc > highLog2fc & negLog10p > negLog10pValue,   1L,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1L, 0L)
    )
  )

  # 5) Build a semantic factor column `category`:
  #    -1 → "down",  0 → "neutral",  1 → "up"
  data$category <- factor(
    data$color_flag,
    levels = c(-1L, 0L, 1L),
    labels = c("down", "neutral", "up"),
    ordered = TRUE
  )

  methods::new("ClearScatterplot", data = data, plot = NULL)
}

#' Internal: run DE per “cell”
#'
#' Performs either limma::voom + lmFit (if count data) or lmFit directly,
#' then eBayes/topTable to get logFC & P‐values for one “cell” (subset of samples).
#'
#' @param expr        numeric matrix of expression values for that cell (features × samples)
#' @param meta        data.frame of metadata rows for those samples
#' @param groupColumn name of the column in `meta` that defines the two groups
#' @param design      design matrix built from `meta[[groupColumn]]`
#' @param dataType    "count" or "continuous"
#' @param cell        list with two elements: SampleType (string) and timePoint (string or NA)
#' @return data.frame with `log2fc`, `negLog10p`, `regulation`, `SampleType`, `timePoint`
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
#' Builds a ClearScatterplot by extracting an assay from a MultiAssayExperiment,
#' filtering low‐variance features, merging in missing top‐level metadata if needed,
#' splitting into “cells” by SampleType (and optional timePoint), running per‐cell DE,
#' and binding all results together.
#'
#' @param mae          a MultiAssayExperiment object
#' @param assayName    string: the assay to extract (must be in `experiments(mae)`)
#' @param groupColumn  string: name of the column in colData (or top‐level colData) defining the groups
#' @param sampleType   string: name of the column in colData (or top‐level colData) used as the facet by X
#' @param timepoint    string or NULL: optional column in colData used as the facet by Y
#' @param dataType     one of c("auto","continuous","count"); "auto" infers based on integer‐ness
#' @param vectorized   one of c("auto","perCell","vectorized"); controls BiocParallel usage
#' @param parallel     logical: if FALSE, disables parallel execution entirely
#' @param BPPARAM      a BiocParallelParam object (e.g., bpparam())
#' @param var_quantile numeric in [0,1]: quantile threshold for filtering low‐variance features
#' @return A `ClearScatterplot` S4 object with a combined DE‐results data.frame
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

  # 1) Ensure the assay exists
  assays <- names(MultiAssayExperiment::experiments(mae))
  if (!assayName %in% assays) {
    stop("Assay '", assayName, "' not found in this MAE.")
  }
  se   <- MultiAssayExperiment::experiments(mae)[[assayName]]
  expr <- SummarizedExperiment::assay(se)

  # 2) Feature‐variance filtering
  rv  <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  # 3) Pull assay‐level metadata (colData)
  meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)

  # 4) If any of groupColumn/sampleType/timepoint are missing, merge from top‐level
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

  # 5) Align expression‐columns and metadata‐rows
  shared <- intersect(colnames(expr), rownames(meta))
  if (length(shared) == 0) {
    stop("No overlapping samples between assay and metadata.")
  }
  expr <- expr[, shared, drop = FALSE]
  meta <- meta[shared, , drop = FALSE]

  # 6) Drop samples with NA in ANY required metadata column
  keep <- rowSums(is.na(meta[needed])) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  # 7) Enforce at least 3 samples per group & at least 6 total
  grp_counts <- table(meta[[groupColumn]])
  if (any(grp_counts < 3)) {
    stop("Each level of '", groupColumn, "' must have ≥ 3 samples.")
  }
  if (ncol(expr) < 6) {
    stop("Too few samples (<6) remain after filtering by NA and variance.")
  }

  # 8) Delegate to the core function, which returns a ClearScatterplot object
  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Constructor: ClearScatterplot_table
#'
#' Build a ClearScatterplot directly from a numeric matrix + matching metadata.
#'
#' @param expr         numeric matrix (features × samples)
#' @param meta         data.frame of sample‐level metadata; row.names must match colnames(expr)
#' @param groupColumn  name of the column in `meta` defining the two groups
#' @param sampleType   name of the column in `meta` used for X‐facet
#' @param timepoint    optional: name of an additional column in `meta` for Y‐facet
#' @param dataType     one of c("auto","continuous","count")—"auto" infers based on integer‐ness
#' @param vectorized   one of c("auto","perCell","vectorized")—controls parallel usage
#' @param parallel     logical: if FALSE, disables parallel entirely
#' @param BPPARAM      a BiocParallelParam object
#' @param var_quantile numeric in [0,1]: quantile threshold for variance filtering
#' @return A `ClearScatterplot` S4 object ready to plot
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

  # 1) Feature‐variance filtering
  rv   <- matrixStats::rowVars(expr, na.rm = TRUE)
  thr  <- stats::quantile(rv, var_quantile, na.rm = TRUE)
  expr <- expr[rv >= thr, , drop = FALSE]

  # 2) Subset metadata to only needed columns, drop NA
  needed <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) needed <- c(needed, timepoint)
  meta2 <- meta[, needed, drop = FALSE]
  keep  <- rowSums(is.na(meta2)) == 0
  expr  <- expr[, keep, drop = FALSE]
  meta2 <- meta2[keep, , drop = FALSE]

  # 3) Enforce at least 3 samples per group & ≥1 feature remains
  grp_counts <- table(meta2[[groupColumn]])
  if (any(grp_counts < 3)) {
    stop("Each level of '", groupColumn, "' must have ≥ 3 samples.")
  }
  if (nrow(expr) < 1) {
    stop("No features remain after variance filtering.")
  }

  # 4) Delegate to core
  .ClearScatterplot_core(
    expr, meta2, groupColumn, sampleType,
    timepoint, dataType,
    if (parallel) vectorized else "perCell",
    BPPARAM
  )
}

#' Internal core for both MAE‐ and matrix‐based constructors
#'
#' Splits into “cells” by (timePoint × SampleType), runs per‐cell DE,
#' binds all results, then calls ClearScatterplot(...) on the combined data.
#'
#' @param expr        numeric matrix (features × samples)
#' @param meta        data.frame; row.names = colnames(expr); must contain grouping/facet info
#' @param groupColumn string giving the grouping column name
#' @param sampleType  string giving the SampleType (X‐facet) column name
#' @param timepoint   string or NULL giving optional Y‐facet column name
#' @param dataType    "auto","continuous","count"
#' @param vectorized  "auto","perCell","vectorized"
#' @param BPPARAM     BiocParallelParam
#' @noRd
.ClearScatterplot_core <- function(
  expr, meta, groupColumn, sampleType,
  timepoint, dataType, vectorized, BPPARAM
) {
  # 1) Drop samples with NA in any of the required columns
  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  keep <- rowSums(is.na(meta[, cols, drop = FALSE])) == 0
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  # 2) Auto‐detect dataType if needed
  if (dataType == "auto") {
    if (all(expr == floor(expr), na.rm = TRUE) && max(expr, na.rm = TRUE) > 30) {
      dataType <- "count"
    } else {
      dataType <- "continuous"
    }
  }

  # 3) Build the “cells” data.frame
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

  # 4) Define per‐cell DE function
  run_cell <- function(i) {
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]
    idx <- if (!is.null(timepoint)) {
      which(meta[[timepoint]] == tp & meta[[sampleType]] == st)
    } else {
      which(meta[[sampleType]] == st)
    }
    # Skip if too few samples or only one group present
    if (length(idx) < 2 || length(unique(meta[[groupColumn]][idx])) < 2) {
      return(NULL)
    }
    ce <- expr[, idx, drop = FALSE]
    cm <- meta[idx, , drop = FALSE]
    # For continuous data, log2‐transform large values
    if (dataType == "continuous" && max(ce, na.rm = TRUE) > 50) {
      ce <- log2(ce + 1)
    }
    design <- stats::model.matrix(~ cm[[groupColumn]])
    if (nrow(cm) <= ncol(design)) return(NULL)
    .run_DE(ce, cm, groupColumn, design, dataType, list(
      SampleType = st, timePoint = tp
    ))
  }

  # 5) Decide on parallelization
  use_par <- (vectorized == "vectorized" ||
              (vectorized == "auto" && nrow(cells) > BiocParallel::bpworkers(BPPARAM)))
  if (use_par) {
    lst <- BiocParallel::bplapply(seq_len(nrow(cells)), run_cell, BPPARAM = BPPARAM)
  } else {
    lst <- lapply(seq_len(nrow(cells)), run_cell)
  }

  # 6) Bind results
  pd <- do.call(rbind, lst)
  if (!nrow(pd)) {
    stop("No DE results to plot (all cells returned zero rows).")
  }
  rownames(pd) <- NULL

  # 7) Build the final ClearScatterplot object (will compute category inside)
  ClearScatterplot(pd)
}

#' Generic: createPlot
#'
#' Build or update the ggplot object inside a `ClearScatterplot` and return it invisibly.
#'
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' Method: createPlot for ClearScatterplot
#'
#' Allows full customization of colors, labels, fonts, etc., with no numeric “magic” hard‐coding.
#'
#' @param object          a ClearScatterplot object
#' @param color1          character color for “up” (category == "up")
#' @param color2          character color for “neutral” (category == "neutral")
#' @param color3          character color for “down” (category == "down")
#' @param xlab            x‐axis label (default: expression(log2~fold~change))
#' @param ylab            y‐axis label (default: expression(-log10~p))
#' @param custom_theme    a ggplot2 theme to add on top (default: NULL)
#' @param point_alpha     numeric 0–1 for point transparency (default: 0.5)
#' @param point_size      numeric for point size (default: 1.75)
#' @param legend_position one of "bottom", "top", "left", "right" (default: "bottom")
#' @param legend_title    title for the color legend (default: NULL)
#' @param legend_labels   character vector of length 3, in order c("down","neutral","up") (default: NULL)
#' @param text_family     font family for all text (default: "sans")
#' @param text_size       base font size for text (default: 10)
#' @param ...             additional arguments passed to ggplot2::labs() (e.g. title, subtitle, caption)
#' @import ggplot2
#' @importFrom dplyr group_by tally
#' @exportMethod createPlot
setMethod("createPlot", "ClearScatterplot", function(object,
  color1          = "cornflowerblue",  # up
  color2          = "grey",            # neutral
  color3          = "indianred",       # down
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

  # Determine facet formula
  has_tp <- "timePoint" %in% names(df) && !all(is.na(df$timePoint))
  facet <- if (has_tp && length(unique(df$timePoint)) > 1) {
    "timePoint ~ SampleType"
  } else {
    ". ~ SampleType"
  }

  # 1) Base ggplot mapping: color by semantic `category`
  p <- ggplot2::ggplot(df, ggplot2::aes(
        x = log2fc,
        y = negLog10p,
        color = category
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

  # 2) Semantic color mapping (no numeric flags):
  #    “up” → color1, “neutral” → color2, “down” → color3
  col_map <- c(
    "up"      = color1,
    "neutral" = color2,
    "down"    = color3
  )
  p <- p + ggplot2::scale_color_manual(
    values = col_map,
    breaks = c("down", "neutral", "up"),  # legend order
    labels = legend_labels
  )

  # 3) Apply user‐supplied theme override (if any)
  if (!is.null(custom_theme)) {
    p <- p + custom_theme
  }

  # 4) Compute and add “up” / “down” count labels
  if (has_tp && length(unique(df$timePoint)) > 1) {
    up <- df[df$category == "up", ] %>%
          dplyr::group_by(timePoint, SampleType) %>%
          dplyr::tally(name = "n1")
    dn <- df[df$category == "down", ] %>%
          dplyr::group_by(timePoint, SampleType) %>%
          dplyr::tally(name = "n2")
  } else {
    up <- df[df$category == "up", ] %>%
          dplyr::group_by(SampleType) %>%
          dplyr::tally(name = "n1")
    dn <- df[df$category == "down", ] %>%
          dplyr::group_by(SampleType) %>%
          dplyr::tally(name = "n2")
  }
  if (nrow(up)) {
    p <- p + ggplot2::geom_text(
      data    = up,
      ggplot2::aes(label = n1),
      x       = Inf,
      y       = Inf,
      hjust   = 1.1,
      vjust   = 1.1,
      color   = color1,            # up‐counts in color1
      family  = text_family,
      size    = text_size / ggplot2::.pt
    )
  }
  if (nrow(dn)) {
    p <- p + ggplot2::geom_text(
      data    = dn,
      ggplot2::aes(label = n2),
      x       = -Inf,
      y       = Inf,
      hjust   = -0.1,
      vjust   = 1.1,
      color   = color3,            # down‐counts in color3
      family  = text_family,
      size    = text_size / ggplot2::.pt
    )
  }

  # 5) Store the finished ggplot object in @plot, return invisibly
  object@plot <- p
  invisible(object)
})

#' Method: show
#'
#' When you type a ClearScatterplot object at the console, automatically
#' build (if needed) and print the ggplot stored in @plot.
#'
#' @param object a ClearScatterplot
#' @rdname ClearScatterplot
#' @exportMethod show
setMethod("show", "ClearScatterplot", function(object) {
  if (is.null(object@plot)) {
    object <- createPlot(object)
  }
  print(object@plot)
  invisible(object)
})
