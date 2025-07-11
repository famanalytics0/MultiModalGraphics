% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ThresholdedScatterplot.R
\name{createPlot}
\alias{createPlot}
\alias{createPlot,ThresholdedScatterplot-method}
\title{Define a Generic Method 'createPlot'}
\usage{
createPlot(object, ...)

\S4method{createPlot}{ThresholdedScatterplot}(
  object,
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
}
\arguments{
\item{object}{A ThresholdedScatterplot object.}

\item{...}{Additional arguments to be passed to the specific method for creating the plot.}

\item{color1}{Color for one category of data points.}

\item{color2}{Color for another category of data points.}

\item{color3}{Color for a third category of data points.}

\item{highLog2fc}{Threshold for high log2 fold change values.}

\item{lowLog2fc}{Threshold for low log2 fold change values.}

\item{negLog10pValue}{Threshold for -log10 p-value.}

\item{expressionDirection}{Direction of gene expression.}

\item{negativeLogPValue}{The name of the column containing the negative log p-values.}

\item{timeVariable}{The variable representing time.}

\item{xAxis}{The x-axis values.}

\item{yAxis}{The y-axis values.}
}
\value{
The result of the specific method for creating the plot, typically a plot object.

The ThresholdedScatterplot object with the plot updated.
}
\description{
This generic function is designed to create plots for objects of various classes.
Specific methods should be defined for different classes to generate the appropriate plots.

Creates a plot based on the ThresholdedScatterplot object data.
}
\examples{
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
}
