# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2"))

#' ClearScatterplot: Class for Enhanced Faceted Volcano Plots
#'
#' S4 class encapsulating data and plotting logic for volcano plots,
#' faceted by sample type and optional timepoint, with count overlays.
#'
#' @slot data data.frame with columns: log2fc, negLog10p, regulation,
#'             SampleType, timePoint (optional), color_flag
#' @slot plot ggplot object for the rendered plot
#' @name ClearScatterplot
#' @docType class
#' @import methods
#' @export
setClass(
  "ClearScatterplot",
  slots = c(
    data = "data.frame",
    plot = "ggplot"
  ),
  validity = function(object) {
    req <- c("log2fc","negLog10p","regulation","SampleType")
    missing <- setdiff(req, colnames(object@data))
    if (length(missing)) {
      return(paste("Missing required columns:", paste(missing, collapse=", ")))
    }
    TRUE
  }
)
    if (length(missing)) {
      return(paste("Missing required columns:", paste(missing, collapse=", ")))
    }
    TRUE
  }
)
)

#' Construct a ClearScatterplot from precomputed data
#'
#' Validates input columns and computes color_flag for visualization.
#'
#' @param data data.frame with at least: log2fc, negLog10p, regulation,
#'        SampleType; optional timePoint column
#' @param highLog2fc numeric threshold for up-regulation (default 0.585)
#' @param lowLog2fc numeric threshold for down-regulation (default -0.585)
#' @param negLog10pValue numeric threshold for significance (default 1.301)
#' @return ClearScatterplot S4 object
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
  highLog2fc = 0.585,
  lowLog2fc = -0.585,
  negLog10pValue = 1.301
) {
  req_cols <- c("log2fc","negLog10p","regulation","SampleType")
  missing <- setdiff(req_cols, colnames(data))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse=", "))
  }
  data$color_flag <- with(data,
    ifelse(log2fc > highLog2fc & negLog10p > negLog10pValue, 1,
    ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1, 0)))
  new("ClearScatterplot", data = data, plot = NULL)
}

#' Build ClearScatterplot from MultiAssayExperiment
#'
#' Runs per-cell limma differential expression and aggregates results
#'
#' @param mae MultiAssayExperiment with assayName
#' @param assayName character name of assay (must exist)
#' @param groupColumn character column in colData for design
#' @param sampleType character column in colData for facet columns
#' @param timepoint character column in colData for facet rows (optional)
#' @return ClearScatterplot object
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom SummarizedExperiment assay
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix complete.cases
#' @export
ClearScatterplot_MAE <- function(
  mae,
  assayName,
  groupColumn = "Group",
  sampleType = "SampleType",
  timepoint = NULL
) {
  if (!inherits(mae, "MultiAssayExperiment")) {
    stop("'mae' must be a MultiAssayExperiment object.")
  }
  assays <- names(MultiAssayExperiment::experiments(mae))
  if (!(assayName %in% assays)) {
    stop("'assayName' must be one of: ", paste(assays, collapse=", "))
  }
  cd <- as.data.frame(MultiAssayExperiment::colData(mae), stringsAsFactors=FALSE)
  req <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) req <- c(req, timepoint)
  missing <- setdiff(req, colnames(cd))
  if (length(missing)) stop("Missing colData cols: ", paste(missing, collapse=", "))
  keep <- rownames(cd)[stats::complete.cases(cd[, req, drop=FALSE])]
  if (length(keep) < 2) stop("Too few samples after dropping NAs.")

  mae2 <- mae[, keep]
  cd2  <- as.data.frame(MultiAssayExperiment::colData(mae2), stringsAsFactors=FALSE)
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
    tp <- cells$timePoint[i]; st <- cells$SampleType[i]
    idx <- if (!is.null(timepoint)) {
      which(cd2[[timepoint]] == tp & cd2[[sampleType]] == st)
    } else {
      which(cd2[[sampleType]] == st)
    }
    if (length(idx) < 2 || length(unique(cd2[[groupColumn]][idx])) < 2) return(NULL)

    expr_c <- expr[, idx, drop=FALSE]
    cd_c   <- cd2[idx, , drop=FALSE]
    design <- stats::model.matrix(~ cd_c[[groupColumn]])
    if (nrow(cd_c) <= ncol(design)) return(NULL)

    fit <- limma::lmFit(expr_c, design) |> limma::eBayes()
    tt  <- limma::topTable(fit, coef=2, number=Inf)

    data.frame(
      log2fc     = tt[['logFC']],
      negLog10p  = -log10(tt[['P.Value']]),
      regulation = ifelse(tt[['logFC']]>0,'up','down'),
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
#' Applies original ggplot styling: geom_point, geom_jitter, theme_bw,
#' and facet_grid with count annotations.
#'
#' @param object ClearScatterplot
#' @param color1 up-regulation color
#' @param color2 neutral color
#' @param color3 down-regulation color
#' @return invisibly returns object with plot slot filled
#' @importFrom ggplot2 ggplot aes geom_point geom_jitter geom_text facet_grid theme_bw theme labs scale_color_manual element_line element_blank
#' @importFrom dplyr group_by tally
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
    facet_formula <- if (!all(is.na(uy)) && length(uy)>1 && length(ux)>1) {
      stats::as.formula(paste(yvar, "~", xvar))
    } else if (!all(is.na(uy)) && length(uy)>1) {
      stats::as.formula(paste(yvar, "~ ."))
    } else {
      stats::as.formula(paste(".~", xvar))
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x=log2fc, y=negLog10p, color=factor(color_flag))) +
      ggplot2::geom_point(alpha=0.5, size=1.75) +
      ggplot2::geom_jitter() +
      ggplot2::labs(x=expression(log2~fold~change), y=expression(-log10~p)) +
      ggplot2::scale_color_manual(values=c(color1,color2,color3)) +
      ggplot2::facet_grid(facet_formula, space="free") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color="grey80"),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill="white", color="black"),
        strip.text       = ggplot2::element_text(size=12,face="bold"),
        axis.title       = ggplot2::element_text(size=12,face="bold"),
        axis.text        = ggplot2::element_text(size=10),
        legend.position  = "bottom"
      )

    cnt_up <- df[df$color_flag==1, ] |> dplyr::group_by(.data[[xvar]], .data[[yvar]]) |> dplyr::tally(name="n1")
    cnt_dn <- df[df$color_flag==-1, ]|> dplyr::group_by(.data[[xvar]], .data[[yvar]]) |> dplyr::tally(name="n2")

    p <- p +
      ggplot2::geom_text(data=cnt_up, ggplot2::aes(label=n1), x=Inf,y=Inf,hjust=1.1,vjust=1.1, color=color1) +
      ggplot2::geom_text(data=cnt_dn, ggplot2::aes(label=n2), x=-Inf,y=Inf,hjust=-0.1,vjust=1.1, color=color3)

    object@plot <- p
    invisible(object)
  }
)

#' Display a ClearScatterplot
#'
#' @param object ClearScatterplot
#' @export
setMethod("show","ClearScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})


