# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2", "negLog10p"))

#' ClearScatterplot: Class for Enhanced Faceted Volcano Scatterplots
#'
#' Provides a unified API for generating volcano plots across sample types
#' and optional timepoints from MultiAssayExperiment objects.
#'
#' @slot data Data frame with columns: log2fc, negLog10p, regulation,
#'             SampleType, timePoint (optional), timeVar, color_flag
#' @slot plot A ggplot object for rendering
#' @name ClearScatterplot
#' @docType class
#' @importFrom methods setClass setGeneric setMethod new
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom SummarizedExperiment assay
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix complete.cases
#' @importFrom ggplot2 ggplot aes geom_point geom_text facet_grid facet_wrap theme scale_color_manual labs
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
#' @param data Data frame of statistics and metadata
#' @param logFoldChange Column name for log2fc
#' @param negativeLogPValue Column name for -log10 p-value
#' @param highLog2fc Threshold for high fold change
#' @param lowLog2fc Threshold for low fold change
#' @param negLog10pValue Significance threshold
#' @param timePointColumn Column storing combined timeVar
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
           timePointColumn = "timeVar")
  {
    stopifnot(is.data.frame(data))
    # Assign flags
    vals <- data[[logFoldChange]]
    pvals <- data[[negativeLogPValue]]
    data[["color_flag"]] <- ifelse(
      vals > highLog2fc & pvals > negLog10pValue, 1,
      ifelse(vals < lowLog2fc & pvals > negLog10pValue, -1, 0)
    )
    # Factorize timeVar if present
    if (timePointColumn %in% names(data)) {
      data[[timePointColumn]] <- as.factor(data[[timePointColumn]])
    }
    new("ClearScatterplot", data = data)
  }
)

#' Build from MultiAssayExperiment
#'
#' @param mae MultiAssayExperiment object
#' @param assayName Assay name to extract
#' @param groupColumn colData column for DE contrast
#' @param sampleType colData column for sample-type facet
#' @param timepoint Optional colData column for timepoint facet
#' @param logFoldChange Name for log2fc output column
#' @param negativeLogPValue Name for negLog10p output column
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
    # Validate assay
    assays <- names(experiments(mae))
    if (is.null(assayName) || !assayName %in% assays) {
      stop("Valid assayName must be one of: ", paste(assays, collapse=", "))
    }
    # Metadata
    cd <- as.data.frame(colData(mae), stringsAsFactors = FALSE)
    req <- c(groupColumn, sampleType)
    if (!all(req %in% colnames(cd))) {
      stop("Missing colData columns: ", paste(setdiff(req, colnames(cd)), collapse=", "))
    }
    # Optional timepoint
    if (!is.null(timepoint) && timepoint %in% colnames(cd)) {
      req <- c(req, timepoint)
      cd <- cd[, req, drop=FALSE]
    } else {
      timepoint <- NULL
    }
    # Drop NA samples
    keep <- rownames(cd)[complete.cases(cd)]
    if (length(keep) < 2) stop("Too few samples after dropping NAs.")
    mae2 <- mae[, keep]
    cd2 <- as.data.frame(colData(mae2), stringsAsFactors = FALSE)
    se  <- experiments(mae2)[[assayName]]
    expr <- assay(se)
    # Build cells: if timepoint, facet grid by both, else only sampleType
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
    }
    # Per-cell DE
    df_list <- lapply(seq_len(nrow(cells)), function(i) {
      if (!is.null(timepoint)) {
        tp <- cells$timePoint[i]
        st <- cells$SampleType[i]
        idx <- which(cd2[[timepoint]]==tp & cd2[[sampleType]]==st)
      } else {
        st <- cells$SampleType[i]
        idx <- which(cd2[[sampleType]]==st)
        tp <- NA_character_
      }
      # require >=2 samples & two groups
      if (length(idx) < 2 || length(unique(cd2[[groupColumn]][idx]))<2) return(NULL)
      expr_c <- expr[, idx, drop=FALSE]
      cd_c   <- cd2[idx, , drop=FALSE]
      design <- model.matrix(~ cd_c[[groupColumn]])
      if (nrow(cd_c) <= ncol(design)) return(NULL)
      fit <- limma::lmFit(expr_c, design) |> limma::eBayes()
      tt  <- limma::topTable(fit, coef=2, number=Inf)
      # Assemble
      df <- data.frame(
        log2fc      = tt[["logFC"]],
        negLog10p   = -log10(tt[["P.Value"]]),
        regulation  = ifelse(tt[["logFC"]]>0, "up", "down"),
        SampleType  = st,
        timePoint   = tp,
        stringsAsFactors = FALSE,
        row.names   = rownames(tt)
      )
      # Combined variable
      if (!is.na(tp)) {
        df[["timeVar"]] <- paste(df[["regulation"]], st, tp, sep="_")
      } else {
        df[["timeVar"]] <- paste(df[["regulation"]], st, sep="_")
      }
      df
    })
    plotdata <- do.call(rbind, df_list)
    rownames(plotdata) <- NULL
    # Return
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
#' @param object ClearScatterplot
#' @param color1, color2, color3 Colors for up/mid/down
#' @param highLog2fc, lowLog2fc Thresholds (currently defaults honored)
#' @param xAxis Facet by SampleType
#' @param yAxis Facet by timeVar
#' @return Updated ClearScatterplot with @plot
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

setMethod(
  "createPlot",
  signature(object="ClearScatterplot"),
  function(object,
           color1="cornflowerblue",
           color2="grey",
           color3="indianred",
           highLog2fc=0.585,
           lowLog2fc=-0.585,
           xAxis="SampleType",
           yAxis="timeVar") {
    df <- object@data
    # Choose facet type
    ux <- unique(df[[xAxis]])
    uy <- unique(df[[yAxis]])
    if (length(uy) > 1 && length(ux) > 1) {
      facet_call <- facet_grid(paste0(yAxis, "~", xAxis), space="free")
    } else if (length(uy) > 1) {
      facet_call <- facet_grid(paste0(yAxis, "~."), space="free")
    } else {
      facet_call <- facet_grid(paste0(".~", xAxis), space="free")
    }
    p <- ggplot(df, aes(x=.data[["log2fc"]], y=.data[["negLog10p"]],
                        color=factor(.data[["color_flag"]]))) +
      geom_point(alpha=0.5) +
      facet_call +
      scale_color_manual(values=c(color1, color2, color3)) +
      theme(strip.text=element_text(face="bold"))
    # Add counts
    up <- df[df$color_flag==1,]
    dn <- df[df$color_flag==-1,]
    cnt_up <- dplyr::group_by(up, .data[[xAxis]], .data[[yAxis]]) |> dplyr::tally(name="n1")
    cnt_dn <- dplyr::group_by(dn, .data[[xAxis]], .data[[yAxis]]) |> dplyr::tally(name="n2")
    p <- p +
      geom_text(data=cnt_up, aes(label=n1), x=Inf, y=Inf, hjust=1.1, vjust=1.1, color=color3) +
      geom_text(data=cnt_dn, aes(label=n2), x=-Inf, y=Inf, hjust=-0.1, vjust=1.1, color=color1) +
      theme(legend.position="bottom")
    object@plot <- p
    object
  }
)

#' Show method for ClearScatterplot
#' @export
setMethod(
  "show",
  signature(object="ClearScatterplot"),
  function(object) {
    if (is.null(object@plot)) object <- createPlot(object)
    print(object@plot)
    invisible(object)
  }
)

