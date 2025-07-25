% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiModalGraphics-package.R
\docType{package}
\name{MultiModalGraphics-package}
\alias{MultiModalGraphics}
\alias{MultiModalGraphics-package}
\title{MultiModalGraphics: An R Package for Graphical Integration of Multimodal Datasets}
\description{
The MultiModalGraphics package provides rich functionalities for combining different data types
into an integrative visual representation. This package facilitates the creation of intuitive visualizations
that merge multiple streams of data into a single coherent graphic.
}
\details{
The MultiModalGraphics package offers the following features:
\itemize{
\item AnnotatedHeatmap: encapsulating a CompexHeatmap object to provide additional functionality
\item CompositeFeatureHeatmap: generating a heatmap grid, where each row presents a composite of information about a group of related features per row, rather than a single feature (unlike typically seen in conventional heatmaps).
\item ThresholdedScatterplot: generating customizable scatter plots that can encompass multiple datasets within a single visualization.
\item generate_data: helper class to provide data to quickly test the 3 classes above.
}
}
\section{Getting Started}{

To get started with the MultiModalGraphics package, you can install it from GitHub:
\preformatted{
# Install the development version from GitHub
# install.packages("devtools")
devtools::install_github("famanalytics0/MultiModalGraphics")
}
}

\examples{
# Example of creating an AnnotatedHeatmap
library(ComplexHeatmap)
library(MultiModalGraphics)
data <- matrix(rnorm(100), ncol = 10)
heatmap <- AnnotatedHeatmap(data, pch_val = 20, unit_val = 2,
                                significant_color = "red",
                                trending_color = "blue",
                                significant_pvalue = 0.05,
                                trending_pvalue = 0.1)
draw(getHeatmapObject(heatmap))

# Example of creating a ThresholdedScatterplot
plotdata <- get_clear_scatterplot_df()
scatterplotObject <- ThresholdedScatterplot(
  data = plotdata,
  logFoldChange = "log2fc",
  timePointColumn = "timePoint",
  timePointLevels = c("T10R1", "T5R1")
)
scattered_plot <- createPlot(
  scatterplotObject,
  color1 = "cornflowerblue",
  color2 = "grey",
  color3 = "indianred",
  highLog2fc = 0.585,
  lowLog2fc = -0.585,
  negLog10pValue = 1.301,
  expressionDirection = "regulation",
  negativeLogPValue = "negLog10p",
  timeVariable = "reg_time_org",
  xAxis = "organ",
  yAxis = "timePoint"
)
print(scattered_plot)

}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/famanalytics0/MultiModalGraphics}
  \item Report bugs at \url{https://github.com/famanalytics0/MultiModalGraphics/issues}
}

}
\author{
\strong{Maintainer}: Seid Muhie \email{yimersm@gmail.com} (\href{https://orcid.org/0000-0003-3409-049X}{ORCID})

}
