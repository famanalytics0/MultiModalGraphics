# Suppress “no visible binding” warnings for variables used in ggplot2
utils::globalVariables(c("Activation_z_score", "neglog10p"))

#’ @title MultifeatureGrid: A Class for 2D Heatmap Visualization
#’ @description
#’ The `MultifeatureGrid` class encapsulates a data‐frame of “signaling” vs. “tissue” 
#’ (or any two categorical axes), along with associated activation z-scores, p-values, 
#’ and item counts.  It provides a single method, `plot_heatmap()`, which builds a fully‐
#’ customized ggplot2 heatmap: tiles are colored by z-score, overplotted points sized by 
#’ item count and colored by –log10(p), and faceting can be applied along any column. 
#’ All computations are fully vectorized, and every input is validated before plotting.
#’ 
#’ @slot data A `data.frame` containing at least the columns:
#’   - `tissue` (factor or character)
#’   - `signaling` (factor or character)
#’   - `Activation_z_score` (numeric)
#’   - a p-value column (e.g. `p`)
#’   - a numeric “count” column (e.g. `number_of_genes`)
#’   - an optional faceting column (default name `timePoint` if present)
#’ @slot title Character; heatmap title.
#’ @slot x_label Character; x‐axis label.
#’ @slot y_label Character; y‐axis label.
#’ @slot logpval_label Character; legend label for –log10(p‐value).
#’ @slot zscore_label Character; legend label for z‐score coloring.
#’ @slot numitems_label Character; legend label for item counts.
#’ @slot color_palette Character; a valid RColorBrewer palette name (length ≥ 3).
#’ @slot breaks Numeric; a strictly increasing vector of breakpoints for the z‐score color scale.
#’ @name MultifeatureGrid
#’ @docType class
#’ @importFrom methods new setClass setGeneric setMethod
#’ @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_gradientn geom_point scale_color_gradient scale_size labs facet_grid theme_bw theme element_text element_rect
#’ @importFrom RColorBrewer brewer.pal
#’ @importFrom grDevices colorRampPalette
#’ @importFrom scales comma_format
#’ @exportClass MultifeatureGrid
#’ @export plot_heatmap
setClass(
  "MultifeatureGrid",
  slots = list(
    data            = "data.frame",
    title           = "character",
    x_label         = "character",
    y_label         = "character",
    logpval_label   = "character",
    zscore_label    = "character",
    numitems_label  = "character",
    color_palette   = "character",
    breaks          = "numeric"
  )
)

#’ @title MultifeatureGrid Constructor
#’ @description
#’ Create a `MultifeatureGrid` object.  Only basic validation is performed here;
#’ most column‐checking happens at plot time.  However, we immediately verify that
#’ `data` is a data.frame and contains the mandatory columns `tissue`, `signaling`, 
#’ and `Activation_z_score`.
#’ 
#’ @param data A `data.frame` with at least the columns:
#’   - `tissue` (factor or character)
#’   - `signaling` (factor or character)
#’   - `Activation_z_score` (numeric)
#’ @param title Character; heatmap title. Default `"Heatmap"`.
#’ @param x_label Character; x‐axis label. Default `"X Label"`.
#’ @param y_label Character; y‐axis label. Default `"Y Label"`.
#’ @param logpval_label Character; legend label for –log10(p‐value). Default `"-log10(p-value)"`.
#’ @param zscore_label Character; legend label for z‐scores. Default `"Activation z-score"`.
#’ @param numitems_label Character; legend label for item counts. Default `"Number of Genes"`.
#’ @param color_palette Character; name of an RColorBrewer palette. Default `"RdYlBu"`.
#’ @param breaks Numeric vector; breakpoints for z‐score color‐mapping. Must be strictly increasing. 
#’   Default `seq(-1,1,0.5)`.
#’   
#’ @return A `MultifeatureGrid` S4 object (no plot yet).
#’ @examples
#’ \dontrun{
#’ df <- data.frame(
#’   tissue = factor(rep(c("T1","T2"), each=4)),
#’   signaling = factor(rep(c("P1","P2","P3","P4"), 2)),
#’   Activation_z_score = runif(8, -2, 2),
#’   p   = runif(8, 0, 0.05),
#’   number_of_genes = sample(1:100, 8),
#’   timePoint = rep(c("Early","Late"), 4),
#’   stringsAsFactors = FALSE
#’ )
#’ mg <- MultifeatureGrid(df, title = "My Heatmap")
#’ }
#’ @export
MultifeatureGrid <- function(
  data,
  title          = "Heatmap",
  x_label        = "X Label",
  y_label        = "Y Label",
  logpval_label  = "-log10(p-value)",
  zscore_label   = "Activation z-score",
  numitems_label = "Number of Genes",
  color_palette  = "RdYlBu",
  breaks         = seq(-1, 1, 0.5)
) {
  ## 1) Basic type‐checking
  if (!inherits(data, "data.frame")) {
    stop("`data` must be a data.frame.")
  }
  required_cols <- c("tissue", "signaling", "Activation_z_score")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("`data` is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  ## 2) Check that Activation_z_score is numeric
  if (!is.numeric(data[["Activation_z_score"]])) {
    stop("Column `Activation_z_score` must be numeric.")
  }
  ## 3) Validate `breaks`
  if (!is.numeric(breaks) || length(breaks) < 2) {
    stop("`breaks` must be a numeric vector of length >= 2.")
  }
  if (any(diff(breaks) <= 0)) {
    stop("`breaks` must be strictly increasing.")
  }
  ## 4) Validate `color_palette` against RColorBrewer
  valid_palettes <- RColorBrewer::brewer.pal.info
  if (!color_palette %in% rownames(valid_palettes)) {
    stop("`color_palette` must be a valid RColorBrewer palette name. ",
         "See `rownames(RColorBrewer::brewer.pal.info)`.")
  }
  ## 5) Everything else must be character of length 1
  for (argnm in c("title","x_label","y_label","logpval_label",
                  "zscore_label","numitems_label","color_palette")) {
    argval <- get(argnm)
    if (!is.character(argval) || length(argval) != 1) {
      stop("`", argnm, "` must be a character string of length 1.")
    }
  }
  ## 6) All checks passed; construct
  methods::new(
    "MultifeatureGrid",
    data            = data,
    title           = title,
    x_label         = x_label,
    y_label         = y_label,
    logpval_label   = logpval_label,
    zscore_label    = zscore_label,
    numitems_label  = numitems_label,
    color_palette   = color_palette,
    breaks          = breaks
  )
}

#’ @title plot_heatmap: Generic for MultifeatureGrid
#’ @description
#’ Generic method for plotting a `MultifeatureGrid`.  Built‐in checks ensure that
#’ the specified p‐value column, count column, and faceting column exist and have
#’ the correct type.  The final ggplot2 call is fully vectorized (no loops).
#’ @param object An object of class `MultifeatureGrid`.
#’ @param ... Additional arguments passed to the specific method.
#’ @export
setGeneric("plot_heatmap", function(object, ...) standardGeneric("plot_heatmap"))

#’ @title plot_heatmap for MultifeatureGrid Objects
#’ @description
#’ Create and display a 2D heatmap with colored tiles (by z‐score), overplotted points
#’ (size = number_of_items, color = –log10(p)), and optional faceting by any column.
#’ This method performs extensive input validation—no loops, all vectorized.
#’
#’ @param object A `MultifeatureGrid` instance.
#’ @param pValueColumn Character(1): name of the column in `object@data` storing raw p‐values. 
#’   Will be transformed to –log10(p).  Default `"p"`.
#’ @param lowColor Character(1): low‐end color for –log10(p) gradient. Default `"yellow"`.
#’ @param highColor Character(1): high‐end color for –log10(p) gradient. Default `"red"`.
#’ @param borderColor Character(1): color of tile borders. Default `"grey60"`.
#’ @param columnForNumber Character(1): name of the column in `data` holding item counts
#’   (e.g. number_of_genes).  These values drive point‐sizes.  Default `"number_of_genes"`.
#’ @param independantVariable Character(1): name of a column to facet columns by.  
#’   If that column is absent, no faceting is applied. Default `"timePoint"`.
#’
#’ @return Invisibly returns the ggplot object, and prints it to the current device.
#’ @export
setMethod("plot_heatmap", signature(object = "MultifeatureGrid"),
          function(object,
                   pValueColumn       = "p",
                   lowColor           = "yellow",
                   highColor          = "red",
                   borderColor        = "grey60",
                   columnForNumber    = "number_of_genes",
                   independantVariable = "timePoint") {
  ## 1) Retrieve slots
  data_df       <- object@data
  title         <- object@title
  xlab          <- object@x_label
  ylab          <- object@y_label
  logp_lbl      <- object@logpval_label
  z_lbl         <- object@zscore_label
  num_lbl       <- object@numitems_label
  palette_name  <- object@color_palette
  breaks_vec    <- object@breaks

  ## 2) Validate data frame columns
  if (!all(c("tissue", "signaling", "Activation_z_score") %in% colnames(data_df))) {
    stop("`data` must contain columns `tissue`, `signaling`, and `Activation_z_score`.")
  }
  if (!is.character(pValueColumn) || length(pValueColumn) != 1) {
    stop("`pValueColumn` must be a single string naming a column in data.")
  }
  if (!(pValueColumn %in% colnames(data_df))) {
    stop("`pValueColumn` '", pValueColumn, "' not found in data.")
  }
  if (!is.numeric(data_df[[pValueColumn]])) {
    stop("Column `", pValueColumn, "` must be numeric (raw p-values).")
  }
  if (any(data_df[[pValueColumn]] <= 0, na.rm = TRUE)) {
    stop("All p-values in `", pValueColumn, "` must be strictly > 0 (no zeros or negatives).")
  }
  if (!is.character(columnForNumber) || length(columnForNumber) != 1) {
    stop("`columnForNumber` must be a single string naming a column in data.")
  }
  if (!(columnForNumber %in% colnames(data_df))) {
    stop("`columnForNumber` '", columnForNumber, "' not found in data.")
  }
  if (!is.numeric(data_df[[columnForNumber]])) {
    stop("Column `", columnForNumber, "` must be numeric (e.g. counts).")
  }
  if (any(data_df[[columnForNumber]] < 0, na.rm = TRUE)) {
    stop("All values in `", columnForNumber, "` must be ≥ 0.")
  }
  if (!is.character(independantVariable) || length(independantVariable) != 1) {
    stop("`independantVariable` must be a single string naming a column in data (or NA).")
  }
  facet_formula <- NULL
  if (independantVariable %in% colnames(data_df)) {
    if (all(is.na(data_df[[independantVariable]]))) {
      warning("`", independantVariable, "` exists but is all NA; facet will be skipped.")
    } else {
      facet_formula <- stats::as.formula(paste(". ~", independantVariable))
    }
  }

  ## 3) Compute –log10(p) and attach
  plot_df <- data_df
  plot_df[["neglog10p"]] <- -log10(data_df[[pValueColumn]])

  ## 4) Build color ramp for Activation_z_score
  brewer_info <- RColorBrewer::brewer.pal.info
  max_colors  <- brewer_info[palette_name, "maxcolors"]
  base_cols   <- RColorBrewer::brewer.pal(max_colors, palette_name)
  color_vec   <- grDevices::colorRampPalette(rev(base_cols))(100)

  ## 5) Build ggplot call
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes_string(
      x    = "tissue",
      y    = "signaling",
      fill = "Activation_z_score"
    )
  ) +
    ggplot2::geom_tile(colour = borderColor) +
    ggplot2::scale_fill_gradientn(
      colours = color_vec,
      breaks  = breaks_vec,
      labels  = scales::comma_format()
    ) +
    ggplot2::geom_point(
      ggplot2::aes_string(
        colour = "neglog10p",
        size   = columnForNumber
      )
    ) +
    ggplot2::scale_color_gradient(
      low  = lowColor,
      high = highColor,
      name = logp_lbl
    ) +
    ggplot2::scale_size(
      range = c(1, 10),
      name  = num_lbl
    ) +
    ggplot2::labs(
      x     = xlab,
      y     = ylab,
      title = title,
      fill  = z_lbl
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey80", size = 0.6),
      axis.text    = ggplot2::element_text(size = 14, face = "bold"),
      axis.title   = ggplot2::element_text(size = 18, face = "bold"),
      plot.title   = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 14, face = "bold", colour = "black")
    )

  ## 6) Add faceting if requested
  if (!is.null(facet_formula)) {
    p <- p + ggplot2::facet_grid(facet_formula, scales = "free_x", space = "free")
  }

  ## 7) Print & return
  print(p)
  invisible(p)
})
