#’ @title ClearScatterplot: Construct a Volcano from a Precomputed Differential Expression Table
#’ @description
#’ Construct a \code{ClearScatterplot} S4 object from a precomputed differential expression table.
#’ The table must contain at least: \code{log2fc}, \code{negLog10p}, \code{regulation}, \code{SampleType}, and optionally \code{timePoint}.
#’ The constructor computes a \code{color_flag} (\code{-1}, \code{0}, \code{1}) and an ordered factor \code{category} (\code{"down"}, \code{"neutral"}, \code{"up"}).
#’ 
#’ @param data A \code{data.frame} with columns:
#’   - \code{log2fc}  (numeric)
#’   - \code{negLog10p} (numeric)
#’   - \code{regulation} (character: \code{"up"} or \code{"down"})
#’   - \code{SampleType}  (character: grouping column)
#’   - (optional) \code{timePoint} (character or factor)
#’ @param highLog2fc Numeric threshold: \code{log2fc > highLog2fc} AND \code{negLog10p > negLog10pValue} ⇒ \code{color_flag = 1} (\code{"up"}).
#’   Default = \code{0.585} (≈ 1.5‐fold).
#’ @param lowLog2fc Numeric threshold: \code{log2fc < lowLog2fc} AND \code{negLog10p > negLog10pValue} ⇒ \code{color_flag = -1} (\code{"down"}).
#’   Default = \code{-0.585}.
#’ @param negLog10pValue Numeric threshold: \code{negLog10p > negLog10pValue} considered significant. Default = \code{1.301} (≈ p < 0.05).
#’ @param dropNA Logical, whether to drop rows with \code{NA} in \code{log2fc} or \code{negLog10p}. Default = \code{TRUE}.
#’ @param ... Additional arguments (reserved for future extensions).
#’ 
#’ @return An S4 object of class \code{ClearScatterplot}, with slots:
#’   - \code{@data}: the processed \code{data.frame} (including \code{color_flag} and \code{category})
#’   - \code{@plot}: \code{NULL} until \code{createPlot()} is called.
#’ 
#’ @examples
#’ \dontrun{
#’ df <- data.frame(
#’   log2fc     = rnorm(200),
#’   negLog10p  = runif(200, 0, 5),
#’   regulation = sample(c("up","down"), 200, TRUE),
#’   SampleType = sample(c("A","B"), 200, TRUE),
#’   stringsAsFactors = FALSE
#’ )
#’ cs <- ClearScatterplot(df)
#’ cs <- createPlot(cs, title = "Volcano Example")
#’ show(cs)
#’ }
#’ @export
ClearScatterplot <- function(
  data,
  highLog2fc     = 0.585,
  lowLog2fc      = -0.585,
  negLog10pValue = 1.301,
  dropNA         = TRUE,
  ...
) {
  # Basic checks
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
  required_cols <- c("log2fc", "negLog10p", "regulation", "SampleType")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  # Check numeric types
  if (!is.numeric(data$log2fc) || !is.numeric(data$negLog10p)) {
    stop("Columns `log2fc` and `negLog10p` must be numeric.")
  }
  # Optionally drop NA rows
  if (dropNA) {
    bad <- is.na(data$log2fc) | is.na(data$negLog10p)
    if (any(bad)) {
      warning(sum(bad), " row(s) removed due to NA in log2fc or negLog10p.")
      data <- data[!bad, , drop = FALSE]
    }
  }
  if (nrow(data) == 0) {
    stop("No rows left after dropping NAs in `log2fc` or `negLog10p`.")
  }
  # Compute color_flag
  data$color_flag <- with(data,
    ifelse(
      log2fc >  highLog2fc & negLog10p > negLog10pValue, 1,
    ifelse(
      log2fc <  lowLog2fc  & negLog10p > negLog10pValue, -1,
      0
    )))
  # Create ordered factor
  data$category <- factor(data$color_flag, levels = c(-1, 0, 1),
                          labels = c("down", "neutral", "up"))
  # Construct S4 object
  methods::new("ClearScatterplot", data = data, plot = NULL)
}


#’ @title createPlot for ClearScatterplot
#’ @description
#’ Build a ggplot2 volcano plot from a \code{ClearScatterplot} object, with faceting by \code{SampleType} and optional \code{timePoint}.
#’ @param object A \code{ClearScatterplot} object whose \code{@data} slot contains the processed DE table.
#’ @param color_up Color for \code{category == "up"}. Default = \code{"indianred"}.
#’ @param color_neutral Color for \code{category == "neutral"}. Default = \code{"grey"}.
#’ @param color_down Color for \code{category == "down"}. Default = \code{"cornflowerblue"}.
#’ @param xlab Label for x-axis. Default = \code{expression(log2~fold~change)}.
#’ @param ylab Label for y-axis. Default = \code{expression(-log10~p)}.
#’ @param legend_position Position of legend: \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}. Default = \code{"bottom"}.
#’ @param legend_title Title for the color legend. Default = \code{NULL}.
#’ @param legend_labels Character vector of length 3 for labels of \code{c("down","neutral","up")}. Default = \code{NULL} (uses factor levels).
#’ @param text_family Font family for all text. Default = \code{"sans"}.
#’ @param text_size Base font size. Default = \code{10}.
#’ @param point_size Size of scatter points. Default = \code{1.75}.
#’ @param point_alpha Transparency of points (0–1). Default = \code{0.5}.
#’ @param custom_theme A ggplot2 theme to apply after \code{theme_bw()}. Default = \code{NULL}.
#’ @param ... Additional arguments passed to \code{ggplot2::labs()}, e.g., \code{title}, \code{subtitle}, \code{caption}.
#’ @return The same \code{ClearScatterplot} object, with its \code{@plot} slot updated to the ggplot2 object.
#’ @export
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

#’ @rdname createPlot
setMethod("createPlot", "ClearScatterplot", function(
  object,
  color_up        = "indianred",
  color_neutral   = "grey",
  color_down      = "cornflowerblue",
  xlab            = expression(log2~fold~change),
  ylab            = expression(-log10~p),
  legend_position = "bottom",
  legend_title    = NULL,
  legend_labels   = NULL,
  text_family     = "sans",
  text_size       = 10,
  point_size      = 1.75,
  point_alpha     = 0.5,
  custom_theme    = NULL,
  ...
) {
  df <- object@data
  # Determine if timePoint exists and is non-NA
  has_tp <- "timePoint" %in% names(df) && !all(is.na(df$timePoint))
  facet  <- if (has_tp && length(unique(df$timePoint)) > 1) {
    "timePoint ~ SampleType"
  } else {
    ". ~ SampleType"
  }
  # Base ggplot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = negLog10p, color = factor(color_flag))) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      color = legend_title,
      ...
    ) +
    ggplot2::scale_color_manual(
      values = c(color_down, color_neutral, color_up),
      labels = legend_labels
    ) +
    ggplot2::facet_grid(stats::as.formula(facet), space = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text             = ggplot2::element_text(family = text_family, size = text_size),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text       = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
      axis.title       = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
      axis.text        = ggplot2::element_text(size = text_size, family = text_family),
      legend.position  = legend_position,
      legend.title     = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
      legend.text      = ggplot2::element_text(size = text_size, family = text_family),
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_blank()
    )
  if (!is.null(custom_theme)) {
    p <- p + custom_theme
  }
  # Counts of up/down per facet
  if (has_tp && length(unique(df$timePoint)) > 1) {
    up_df <- df[df$color_flag == 1, ] %>%
      dplyr::group_by(timePoint, SampleType) %>%
      dplyr::tally(name = "n_up")
    dn_df <- df[df$color_flag == -1, ] %>%
      dplyr::group_by(timePoint, SampleType) %>%
      dplyr::tally(name = "n_down")
  } else {
    up_df <- df[df$color_flag == 1, ] %>%
      dplyr::group_by(SampleType) %>%
      dplyr::tally(name = "n_up")
    dn_df <- df[df$color_flag == -1, ] %>%
      dplyr::group_by(SampleType) %>%
      dplyr::tally(name = "n_down")
  }
  if (nrow(up_df) > 0) {
    p <- p + ggplot2::geom_text(
      data = up_df,
      ggplot2::aes(label = n_up),
      x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
      color = color_up, family = text_family,
      size = text_size / ggplot2::.pt
    )
  }
  if (nrow(dn_df) > 0) {
    p <- p + ggplot2::geom_text(
      data = dn_df,
      ggplot2::aes(label = n_down),
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
      color = color_down, family = text_family,
      size = text_size / ggplot2::.pt
    )
  }
  object@plot <- p
  invisible(object)
})

#’ @title Show method for ClearScatterplot
#’ @param object A \code{ClearScatterplot} object with \code{@plot} slot built.
#’ @rdname createPlot
#’ @export
setMethod("show", "ClearScatterplot", function(object) {
  if (is.null(object@plot)) {
    stop("Plot has not been created yet. Run `createPlot(object)` first.")
  }
  print(object@plot)
  invisible(NULL)
})
