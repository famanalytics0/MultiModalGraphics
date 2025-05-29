# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2", "negLog10p"))

#' ClearScatterplot: Class for Enhanced Faceted Volcano Scatterplots
#'
#' Provides a unified API for generating volcano plots across sample types
#' and optional timepoints from MultiAssayExperiment objects.
#'
#' @slot data A data frame with columns: log2fc, negLog10p, regulation,
#'             SampleType, timePoint (optional), timeVar, color_flag
#' @slot plot A ggplot object storing the rendered plot
#' @name ClearScatterplot
#' @docType class
#' @importFrom methods setClass setGeneric setMethod new
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom SummarizedExperiment assay
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix complete.cases
#' @importFrom ggplot2 ggplot aes geom_point geom_text facet_grid theme scale_color_manual labs
#' @importFrom dplyr group_by tally
#' @export
setClass(
  "ClearScatterplot",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  )
)

#' Constructor for ClearScatterplot
#'
#' @param data A data.frame containing computed statistics and metadata
#' @param logFoldChange Column name for log2 fold change values
#' @param negativeLogPValue Column name for -log10 p-values
#' @param highLog2fc Threshold above which fold change is considered high
#' @param lowLog2fc Threshold below which fold change is considered low
#' @param negLog10pValue Threshold above which p-value is significant
#' @param timePointColumn Column name for the combined time variable
#' @return A ClearScatterplot object
#' @export
setGeneric(
  "ClearScatterplot",
  function(data,
           logFoldChange = "log2fc",
           negativeLogPValue = "negLog10p",
           highLog2fc = 0.585,
           lowLog2fc = -0.585,
           negLog10pValue = 1.301,
           timePointColumn = "timePoint") {
    stopifnot(is.data.frame(data))
    # Compute color_flag
    vals <- data[[logFoldChange]]
    pvals <- data[[negativeLogPValue]]
    data[["color_flag"]] <- ifelse(
      vals > highLog2fc & pvals > negLog10pValue, 1,
      ifelse(vals < lowLog2fc & pvals > negLog10pValue, -1, 0)
    )
    # Factorize timepoint if present
    if (timePointColumn %in% names(data)) {
      data[[timePointColumn]] <- as.factor(data[[timePointColumn]])
    }
    new("ClearScatterplot", data = data)
  }
)

#' Build ClearScatterplot from MultiAssayExperiment
#'
#' @param mae A MultiAssayExperiment object
#' @param assayName Name of the assay (SummarizedExperiment) to use
#' @param groupColumn Column in colData defining groups for DE
#' @param sampleType Column in colData for sample-type faceting
#' @param timepoint Optional column in colData for timepoint faceting
#' @param logFoldChange Name for output log2fc column
#' @param negativeLogPValue Name for output negLog10p column
#' @return A ClearScatterplot object
#' @export
setGeneric(
  "ClearScatterplot_MAE",
  function(mae,
           assayName = NULL,
           groupColumn = "Group",
           sampleType = "SampleType",
           timepoint = NULL,
           logFoldChange = "log2fc",
           negativeLogPValue = "negLog10p")
    standardGeneric("ClearScatterplot_MAE")
)

setMethod(
  "ClearScatterplot_MAE",
  signature(mae = "MultiAssayExperiment"),
  function(mae,
           assayName = NULL,
           groupColumn = "Group",
           sampleType = "SampleType",
           timepoint = NULL,
           logFoldChange = "log2fc",
           negativeLogPValue = "negLog10p") {
    # Validate assayName
    assays <- names(experiments(mae))
    if (is.null(assayName) || !assayName %in% assays) {
      stop("assayName must be one of: ", paste(assays, collapse=", "))
    }
    # Extract colData
    cd <- as.data.frame(colData(mae), stringsAsFactors = FALSE)
    # Required columns
    req <- c(groupColumn, sampleType)
    if (!all(req %in% colnames(cd))) {
      stop("Missing colData columns: ", paste(setdiff(req, colnames(cd)), collapse=", "))
    }
    # Include timepoint if provided
    if (!is.null(timepoint) && timepoint %in% colnames(cd)) {
      req <- c(req, timepoint)
    } else {
      timepoint <- NULL
    }
    # Subset samples with complete metadata
    keep <- rownames(cd)[complete.cases(cd[, req, drop = FALSE])]
    if (length(keep) < 2) stop("Too few samples after dropping NAs.")
    mae2 <- mae[, keep]
    cd2  <- as.data.frame(colData(mae2), stringsAsFactors = FALSE)
    se   <- experiments(mae2)[[assayName]]
    expr <- assay(se)
    # Build combination grid
    if (!is.null(timepoint)) {
      cells <- expand.grid(
        timePoint  = unique(cd2[[timepoint]]),
        SampleType = unique(cd2[[sampleType]]),
        stringsAsFactors = FALSE
      )
    } else {
      cells <- data.frame(
        SampleType = unique(cd2[[sampleType]]),
        stringsAsFactors = FALSE
      )
      cells$timePoint <- NA_character_
    }
    # Run per-cell DE
    df_list <- lapply(seq_len(nrow(cells)), function(i) {
      tp <- cells$timePoint[i]
      st <- cells$SampleType[i]
      if (!is.null(timepoint)) {
        idx <- which(cd2[[timepoint]] == tp & cd2[[sampleType]] == st)
      } else {
        idx <- which(cd2[[sampleType]] == st)
      }
      # Minimum samples & groups
      if (length(idx) < 2 || length(unique(cd2[[groupColumn]][idx])) < 2) return(NULL)
      expr_c <- expr[, idx, drop = FALSE]
      cd_c   <- cd2[idx, , drop = FALSE]
      design <- model.matrix(~ cd_c[[groupColumn]])
      if (nrow(cd_c) <= ncol(design)) return(NULL)
      fit <- limma::lmFit(expr_c, design) |> limma::eBayes()
      tt  <- limma::topTable(fit, coef = 2, number = Inf)
      # Assemble results
      df <- data.frame(
        log2fc     = tt[["logFC"]],
        negLog10p  = -log10(tt[["P.Value"]]),
        regulation = ifelse(tt[["logFC"]] > 0, "up", "down"),
        SampleType = st,
        timePoint  = tp,
        stringsAsFactors = FALSE,
        row.names  = rownames(tt)
      )
      # Combined time variable
      df[["timeVar"]] <- if (!is.na(tp)) {
        paste(df[["regulation"]], st, tp, sep = "_")
      } else {
        paste(df[["regulation"]], st, sep = "_")
      }
      df
    })
    plotdata <- do.call(rbind, df_list)
    rownames(plotdata) <- NULL
    # Construct object
    ClearScatterplot(
      data = plotdata,
      logFoldChange = logFoldChange,
      negativeLogPValue = negativeLogPValue,
      timePointColumn = "timeVar"
    )
  }
)

#' Generate Plot from ClearScatterplot
#'
#' @param object A ClearScatterplot object
#' @param color1 Color for high up-regulation
#' @param color2 Color for neutral
#' @param color3 Color for high down-regulation
#' @param highLog2fc Threshold for high fold-change
#' @param lowLog2fc Threshold for low fold-change
#' @param xAxis Column for facet columns (sampleType)
#' @param yAxis Column for facet rows (timePoint)
#' @return The same ClearScatterplot object with @plot populated
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

setMethod(
  "createPlot",
  signature(object = "ClearScatterplot"),
  function(object,
           color1 = "cornflowerblue",
           color2 = "grey",
           color3 = "indianred",
           highLog2fc = 0.585,
           lowLog2fc = -0.585,
           xAxis = "SampleType",
           yAxis = "timePoint") {
    df <- object@data
    ux <- unique(df[[xAxis]])
    uy <- unique(df[[yAxis]])
    # Dynamic faceting
    facet_formula <- if (length(uy) > 1 && length(ux) > 1) {
      stats::as.formula(paste(yAxis, "~", xAxis))
    } else if (length(uy) > 1) {
      stats::as.formula(paste(yAxis, "~ ."))
    } else {
      stats::as.formula(paste(".~", xAxis))
    }
    # Base plot
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data[["log2fc"]],
      y = .data[["negLog10p"]],
      color = factor(.data[["color_flag"]])
    )) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::facet_grid(facet_formula, space = "free") +
      ggplot2::scale_color_manual(values = c(color1, color2, color3)) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
    # Overlay counts
    df_up <- subset(df, color_flag == 1)
    df_dn <- subset(df, color_flag == -1)
    cnt_up <- df_up |>
      dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) |> dplyr::tally(name = "n1")
    cnt_dn <- df_dn |>
      dplyr::group_by(.data[[xAxis]], .data[[yAxis]]) |> dplyr::tally(name = "n2")
    p <- p +
      ggplot2::geom_text(
        data = cnt_up,
        ggplot2::aes(label = n1),
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 1.1,
        color = color1
      ) +
      ggplot2::geom_text(
        data = cnt_dn,
        ggplot2::aes(label = n2),
        x = -Inf, y = Inf,
        hjust = -0.1, vjust = 1.1,
        color = color3
      )
    object@plot <- p
    invisible(object)
  }
)

#' Show method for ClearScatterplot
#' @export
setMethod(
  "show",
  signature(object = "ClearScatterplot"),
  function(object) {
    if (is.null(object@plot)) object <- createPlot(object)
    print(object@plot)
    invisible(object)
  }
)


