% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CompositeFeatureHeatmap.R
\name{plot_heatmap,CompositeFeatureHeatmap-method}
\alias{plot_heatmap,CompositeFeatureHeatmap-method}
\title{Plot Heatmap for CompositeFeatureHeatmap Objects}
\usage{
\S4method{plot_heatmap}{CompositeFeatureHeatmap}(
  object,
  pValueColumn = "p",
  lowColor = "yellow",
  highColor = "red",
  borderColor = "grey60",
  columnForNumber = "number_of_genes",
  independantVariable = "timePoint"
)
}
\arguments{
\item{object}{An object of class \code{CompositeFeatureHeatmap}. This object contains
the data and configuration for the heatmap to be plotted.}

\item{pValueColumn}{A character string specifying the name of the column in the
data frame that contains p-values. These values are used to calculate
and plot the negative log10 of the p-values. Default is "p".}

\item{lowColor}{A character string specifying the color to use for the low
end of the color gradient for p-values. Default is "yellow".}

\item{highColor}{A character string specifying the color to use for the high
end of the color gradient for p-values. Default is "red".}

\item{borderColor}{A character string specifying the color to use for the border
of each tile in the heatmap. Default is "grey60".}

\item{columnForNumber}{A character string specifying the name of the column in the
data frame that contains the number of items (e.g., genes) for each feature.
This value is used to determine the size of the points plotted on the heatmap.
Default is "number_of_genes".}

\item{independantVariable}{A character string specifying the name of the column used
as an independent variable for faceting the plot into multiple panels, one for each
level of the independent variable. Default is "timePoint".}
}
\value{
A ggplot object representing the heatmap, which is also printed to the current
plotting device.
}
\description{
This method generates and plots a heatmap for \code{CompositeFeatureHeatmap} objects. It allows
customization of the heatmap's appearance based on several parameters.
}
\examples{
data <- data.frame(tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
                   signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
                   Activation_z_score = runif(8, -2, 2),
                   p = runif(8, 0, 0.05),
                   number_of_genes = sample(1:100, 8),
                   timePoint = rep(c("Time1", "Time2"), 4))
mg <- CompositeFeatureHeatmap(data)
plot_heatmap(mg)
}
