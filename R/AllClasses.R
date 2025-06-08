#' S4 Class Definitions for MultiModalGraphics
#'
#' This file defines all core S4 classes used by the MultiModalGraphics package.
#' Keep all S4 setClass definitions here to ensure reproducibility, extensibility,
#' and correct method dispatch across the package.
#'
#' @section Core classes:
#' \describe{
#'   \item{ClearScatterplot}{Volcano scatterplots for precomputed or computed DE tables}
#'   \item{InformativeHeatmap}{Layered/annotated heatmaps from feature matrices}
#'   \item{MultifeatureGrid}{Multi-panel feature-by-feature grid visualizations}
#' }
#' @name AllClasses
#' @docType package
NULL

# ==== ClearScatterplot S4 Class ====

#' @title ClearScatterplot S4 Class
#' @description
#' Stores a precomputed DE result (or output from a MultiAssayExperiment/matrix) and its ggplot2 visualization.
#' @slot data data.frame. The processed DE table.
#' @slot plot ANY. The plot object (ggplot2), constructed on demand.
#' @export
setClass(
  Class = "ClearScatterplot",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  ),
  prototype = list(
    data = data.frame(),
    plot = NULL
  )
)

# ==== InformativeHeatmap S4 Class ====

#' @title InformativeHeatmap S4 Class
#' @description
#' Stores the data and ComplexHeatmap object for an annotated, layered heatmap.
#' @slot data data.frame. Matrix or feature metadata used in the heatmap.
#' @slot plot ANY. The ComplexHeatmap/Heatmap object, constructed on demand.
#' @export
setClass(
  Class = "InformativeHeatmap",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  ),
  prototype = list(
    data = data.frame(),
    plot = NULL
  )
)

# ==== MultifeatureGrid S4 Class ====

#' @title MultifeatureGrid S4 Class
#' @description
#' Stores the data and plot object for a grid of per-feature or per-module plots.
#' @slot data data.frame. Data for grid construction.
#' @slot plot ANY. ggplot2 patchwork or list of ggplots.
#' @export
setClass(
  Class = "MultifeatureGrid",
  slots = c(
    data = "data.frame",
    plot = "ANY"
  ),
  prototype = list(
    data = data.frame(),
    plot = NULL
  )
)
