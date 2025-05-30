# Suppress warnings for global variables
utils::globalVariables(c(".data", "n1", "n2"))

#' S4 Class: ClearScatterplot
#'
#' Encapsulates data and plotting logic for volcano plots
#' with flexible faceting.
#'
#' @slot data   data.frame with plot data
#' @slot plot   ggplot object
#' @exportClass ClearScatterplot
setClass(
  "ClearScatterplot",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  ),
  validity = function(object) {
    req <- c("log2fc", "negLog10p", "regulation", "SampleType")
    miss <- setdiff(req, names(object@data))
    if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
    TRUE
  }
)

#' Generic: ClearScatterplot
#' @export
setGeneric(
  "ClearScatterplot",
  function(data, ...) standardGeneric("ClearScatterplot")
)

#' Constructor: ClearScatterplot
#'
#' @param data            data.frame
#' @param highLog2fc      threshold (default 0.585)
#' @param lowLog2fc       threshold (default -0.585)
#' @param negLog10pValue  threshold (default 1.301)
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
  if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
  data$color_flag <- with(data,
    ifelse(log2fc > highLog2fc & negLog10p > negLog10pValue, 1,
      ifelse(log2fc < lowLog2fc & negLog10p > negLog10pValue, -1, 0)
    )
  )
  new("ClearScatterplot", data = data, plot = NULL)
}

# Internal: run DE for one cell
.run_DE <- function(expr, meta, groupColumn, designMat, dataType, cell) {
  fit <- if (dataType == "count") {
    voom <- limma::voom(expr, designMat, plot = FALSE)
    limma::lmFit(voom, designMat)
  } else {
    limma::lmFit(expr, designMat)
  }
  fit <- limma::eBayes(fit)
  tt  <- limma::topTable(fit, coef=2, number=Inf)
  if (!nrow(tt)) return(NULL)
  data.frame(
    log2fc     = tt$logFC,
    negLog10p  = -log10(tt$P.Value),
    regulation = ifelse(tt$logFC>0, "up", "down"),
    SampleType = cell$SampleType,
    timePoint  = cell$timePoint,
    row.names  = rownames(tt),
    stringsAsFactors = FALSE
  )
}

#' Constructor: from MultiAssayExperiment
#' @export
ClearScatterplot_MAE <- function(
  mae, assayName,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto","continuous","count"),
  vectorized  = c("auto","perCell","vectorized"),
  BPPARAM     = BiocParallel::bpparam(),
  var_quantile=0.75
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)

  full <- SummarizedExperiment::assay(
    MultiAssayExperiment::experiments(mae)[[assayName]]
  )
  rv   <- matrixStats::rowVars(full, na.rm=TRUE)
  feat <- which(rv >= quantile(rv, var_quantile, na.rm=TRUE))
  mae  <- subsetByRow(mae, IRanges::IntegerList(setNames(list(feat), assayName)))

  meta <- as.data.frame(
    SummarizedExperiment::colData(mae),
    stringsAsFactors = FALSE
  )
  expr <- SummarizedExperiment::assay(
    MultiAssayExperiment::experiments(mae)[[assayName]]
  )

  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  md   <- meta[, cols, drop=FALSE]
  keep <- rownames(md)[apply(md,1, function(r) all(!is.na(r)))]
  meta <- meta[keep,,drop=FALSE]
  expr <- expr[,keep,drop=FALSE]

  tab <- table(meta[[groupColumn]])
  if (any(tab<3)) stop("<3 samples per group")
  if (length(keep)<6 ) stop("<6 samples total")

  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType, vectorized, BPPARAM
  )
}

#' Constructor: from matrix + metadata
#' @export
ClearScatterplot_table <- function(
  expr, meta,
  groupColumn = "Group",
  sampleType  = "SampleType",
  timepoint   = NULL,
  dataType    = c("auto","continuous","count"),
  vectorized  = c("auto","perCell","vectorized"),
  BPPARAM     = BiocParallel::bpparam(),
  var_quantile=0.75
) {
  dataType   <- match.arg(dataType)
  vectorized <- match.arg(vectorized)
  stopifnot(is.matrix(expr), is.data.frame(meta))

  rv   <- matrixStats::rowVars(expr, na.rm=TRUE)
  feat <- which(rv >= quantile(rv,var_quantile,na.rm=TRUE))
  expr <- expr[feat,,drop=FALSE]

  cols <- c(groupColumn, sampleType)
  if (!is.null(timepoint)) cols <- c(cols, timepoint)
  md   <- meta[, cols, drop=FALSE]
  keep <- rownames(md)[apply(md,1, function(r) all(!is.na(r)))]
  meta <- meta[keep,,drop=FALSE]
  expr <- expr[,keep,drop=FALSE]

  tab <- table(meta[[groupColumn]])
  if (any(tab<3)) stop("<3 samples per group")
  if (length(keep)<6 ) stop("<6 samples total")

  .ClearScatterplot_core(
    expr, meta, groupColumn, sampleType,
    timepoint, dataType, vectorized, BPPARAM
  )
}

# Internal core shared by both
.ClearScatterplot_core <- function(
  expr, meta, groupColumn, sampleType,
  timepoint, dataType, vectorized, BPPARAM
) {
  keep <- rownames(meta)[complete.cases(meta[,c(groupColumn,sampleType,timepoint),drop=FALSE])]
  expr <- expr[,keep,drop=FALSE]
  meta <- meta[keep,,drop=FALSE]

  if (dataType=="auto") {
    dataType <- if(all(expr==floor(expr))&&max(expr,na.rm=TRUE)>30) "count" else "continuous"
  }

  cells <- if (!is.null(timepoint)) {
    expand.grid(
      timePoint=unique(meta[[timepoint]]),
      SampleType=unique(meta[[sampleType]]),
      stringsAsFactors=FALSE
    )
  } else {
    data.frame(
      SampleType=unique(meta[[sampleType]]),
      timePoint=NA_character_,
      stringsAsFactors=FALSE
    )
  }

  run_cell <- function(i) {
    tp <- cells$timePoint[i]
    st <- cells$SampleType[i]
    idx<- if(!is.null(timepoint)) {
      which(meta[[timepoint]]==tp & meta[[sampleType]]==st)
    } else {
      which(meta[[sampleType]]==st)
    }
    if (length(idx)<2||length(unique(meta[[groupColumn]][idx]))<2) return(NULL)
    .run_DE(
      expr[,idx,drop=FALSE],
      meta[idx,,drop=FALSE],
      groupColumn,
      model.matrix(~meta[idx,groupColumn]),
      dataType,
      list(SampleType=st,timePoint=tp)
    )
  }

  df_list <- if(vectorized=="vectorized"||
                (vectorized=="auto"&&nrow(cells)>BiocParallel::bpworkers(BPPARAM))) {
    BiocParallel::bplapply(seq_len(nrow(cells)),run_cell,BPPARAM=BPPARAM)
  } else {
    lapply(seq_len(nrow(cells)),run_cell)
  }

  plotdata<-do.call(rbind,df_list)
  if(!nrow(plotdata)) stop("No DE results to plot")
  rownames(plotdata)<-NULL
  ClearScatterplot(plotdata)
}

#' Generic: createPlot
#' @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#' Method: createPlot for ClearScatterplot
#' @exportMethod createPlot
setMethod(
  "createPlot",
  "ClearScatterplot",
  function(object, color1="cornflowerblue", color2="grey", color3="indianred") {
    df <- object@data
    x  <- "SampleType"; y <- "timePoint"
    ux <- unique(df[[x]]); uy <- unique(df[[y]])
    facet <- if(length(ux)>1&&length(uy)>1) {
      stats::as.formula(paste(y,"~",x))
    } else if(length(uy)>1) {
      stats::as.formula(paste(y,"~ ."))
    } else {
      stats::as.formula(paste(".~",x))
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x=log2fc,y=negLog10p,color=factor(color_flag))) +
      ggplot2::geom_point(alpha=0.5,size=1.75) +
      ggplot2::geom_jitter() +
      ggplot2::labs(x=expression(log2~fold~change),y=expression(-log10~p)) +
      ggplot2::scale_color_manual(values=c(color1,color2,color3)) +
      ggplot2::facet_grid(facet, space="free") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major=ggplot2::element_line(color="grey80"),
        panel.grid.minor=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="white",color="black"),
        strip.text=ggplot2::element_text(size=12,face="bold"),
        axis.title=ggplot2::element_text(size=12,face="bold"),
        axis.text=ggplot2::element_text(size=10),
        legend.position="bottom"
      )

    cnt_up <- df[df$color_flag==1,] |>
      dplyr::group_by(.data[[x]],.data[[y]]) |> dplyr::tally(name="n1")
    cnt_dn <- df[df$color_flag==-1,]|>
      dplyr::group_by(.data[[x]],.data[[y]]) |> dplyr::tally(name="n2")

    p <- p +
      ggplot2::geom_text(
        data=cnt_up, ggplot2::aes(label=n1),
        x=Inf,y=Inf,hjust=1.1,vjust=1.1,color=color1
      ) +
      ggplot2::geom_text(
        data=cnt_dn, ggplot2::aes(label=n2),
        x=-Inf,y=Inf,hjust=-0.1,vjust=1.1,color=color3
      )

    object@plot <- p
    invisible(object)
  }
)

#' Method: show for ClearScatterplot
#' @exportMethod show
setMethod("show","ClearScatterplot", function(object) {
  if (is.null(object@plot)) object <- createPlot(object)
  print(object@plot)
  invisible(object)
})






