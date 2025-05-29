# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2"))

#' ClearScatterplot: Class for Enhanced Faceted Volcano Plots
#'
#' S4 class encapsulating data and plotting logic for volcano plots,
#' faceted by sample type and optional timepoint, with count overlays.
#'
#' @slot data data.frame with columns: log2fc, negLog10p, regulation,
#'             SampleType, timePoint (optional), color_flag
#' @slot plot  ANY storing the generated ggplot object
#' @name ClearScatterplot
#' @docType class
#' @import methods
#' @importFrom ggplot2 ggplot aes geom_point geom_jitter geom_text facet_grid theme_bw theme labs scale_color_manual element_line element_blank
#' @importFrom dplyr group_by tally
#' @export
setClass(
  "ClearScatterplot",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  ),
  validity = function(object) {
    req_cols <- c("log2fc", "negLog10p", "regulation", "SampleType")
    missing  <- setdiff(req_cols, names(object@data))
    if (length(missing)) {
      stop("Missing required columns: ", paste(missing, collapse = ", "))
    }
    TRUE
  }
)

#' Constructor: ClearScatterplot
#'
#' Validates input data and computes a color flag for plotting thresholds.
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
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  data$color_flag <- with(
    data,
    ifelse(
      log2fc > highLog2fc & negLog10p > negLog10pValue, 1,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1, 0)
    )
  )

  new("ClearScatterplot", data = data, plot = NULL)
}

#' Constructor: ClearScatterplot from MAE
#'
#' Builds ClearScatterplot by running per-cell limma differential
#' expression on a MultiAssayExperiment.
#'
#' @param mae         a MultiAssayExperiment
#' @param assayName   name of the assay to use
#' @param groupColumn colData column for group design
#' @param sampleType  colData column for facet columns
#' @param timepoint   optional colData column for facet rows
#' @return ClearScatterplot instance
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom SummarizedExperiment assay
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix complete.cases
#' @export
ClearScatterplot_MAE <- function(
  mae,
  assayName,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL
) {
  stopifnot(inherits(mae, "MultiAssayExperiment"))
  assays <- names(MultiAssayExperiment::experiments(mae))
  if (!assayName %in% assays) {
    stop("assayName must be one of: ", paste(assays, collapse = ", "))
  }

  cd <- as.data.frame(
    MultiAssayExperiment::colData(mae),
    stringsAsFactors = FALSE
  )
  req <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) req <- c(req, timepoint)
  missing_cd <- setdiff(req, names(cd))
  if (length(missing_cd)) {
    stop("Missing colData columns: ", paste(missing_cd, collapse = ", "))
  }

  keep <- rownames(cd)[complete.cases(cd[req])]
  if (length(keep) < 2) {
    stop("Too few samples after dropping NAs in: ", paste(req, collapse = ", "))
  }

  mae2 <- mae[, keep]
  cd2  <- as.data.frame(
    MultiAssayExperiment::colData(mae2),
    stringsAsFactors = FALSE
  )
  se   <- MultiAssayExperiment::experiments(mae2)[[assayName]]
  expr <- SummarizedExperiment::assay(se)

  cells <- if (!is.null(timepoint)) {
    expand.grid(
      timePoint  = unique(cd2[[timepoint]]),
      SampleType = unique(cd2[[sampleType]]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      SampleType = unique(cd2[[sampleType]]),
      timePoint  = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  df_list <- lapply(seq_len(nrow(cells)), function(i) {
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]

    idx <- if (!is.null(timepoint)) {
      which(
        cd2[[timepoint]]  == tp &
        cd2[[sampleType]] == st
      )
    } else {
      which(cd2[[sampleType]] == st)
    }

    if (length(idx) < 2 || length(unique(cd2[[groupColumn]][idx])) < 2) {
      return(NULL)
    }

    expr_c <- expr[, idx, drop = FALSE]
    cd_c   <- cd2[idx, , drop = FALSE]
    design <- stats::model.matrix(~ cd_c[[groupColumn]])
    if (nrow(cd_c) <= ncol(design)) return(NULL)

    fit <- limma::lmFit(expr_c, design) |> limma::eBayes()
    tt  <- limma::topTable(fit, coef = 2, number = Inf)

    data.frame(
      log2fc     = tt[["logFC"]],
      negLog10p  = -log10(tt[["P.Value"]]),
      regulation = ifelse(tt[["logFC"]] > 0, "up", "down"),
      SampleType = st,
      timePoint  = tp,
      stringsAsFactors = FALSE,
      row.names  = rownames(tt)
    )
  })

  plotdata <- do.call(rbind, df_list)
  rownames(plotdata) <- NULL

  ClearScatterplot(plotdata)
}

#' Render a ClearScatterplot
#'
#' Applies ggplot styling: points, jitter, theme_bw,
#' and facet_grid with count annotations.
#'
#' @param object ClearScatterplot
#' @param color1 up-regulation color
#' @param color2 neutral color
#' @param color3 down-regulation color
#' @return object with plot slot filled
#' @export
setGeneric(
  "createPlot",
  function(object, ...) standardGeneric("createPlot")
)
setMethod(
  "createPlot",
  signature(object = "ClearScatterplot"),
  function(
    object,
    color1 = "cornflowerblue",
    color2 = "grey",
    color3 = "indianred"
  ) {
    df   <- object@data
    xvar <- "SampleType"
    yvar <- "timePoint"
    ux   <- unique(df[[xvar]])
    uy   <- unique(df[[yvar]])
    facet_formula <- if (
      length(ux) > 1 && length(uy) > 1
    ) {
      stats::as.formula(paste(yvar, "~", xvar))
    } else if (length(uy) > 1) {
      stats::as.formula(paste(yvar, "~ ."))
    } else {
      stats::as.formula(paste(".~", xvar))
    }

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x     = log2fc,
        y     = negLog10p,
        color = factor(color_flag)
      )
    ) +
      ggplot2::geom_point(alpha = 0.5, size = 1.75) +
      ggplot2::geom_jitter() +
      ggplot2::labs(
        x = expression(log2 ~ fold ~ change),
        y = expression(-log10 ~ p)
      ) +
      ggplot2::scale_color_manual(
        values = c(color1, color2, color3)
      ) +
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

    cnt_up <- df[df$color_flag == 1, ] |> dplyr::group_by(.data[[xvar]], .data[[yvar]]) |> dplyr::tally(name = "n1")
    cnt_dn <- df[df$color_flag == -1, ] |> dplyr::group_by(.data[[xvar]], .data[[yvar]]) |> dplyr::tally(name = "n2")

    p <- p +
      ggplot2::geom_text(data = cnt_up, ggplot2::aes(label = n1), x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, color = color1) +
      ggplot2::geom_text(data = cnt_dn, ggplot2::aes(label = n2), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, color = color3)

    object@plot <- p
    invisible(object)
  }
)

#' Show method
#' @export
setMethod(
  "show",
  "ClearScatterplot",
  function(object) {
    if (is.null(object@plot)) object <- createPlot(object)
    print(object@plot)
    invisible(object)
  }
)

