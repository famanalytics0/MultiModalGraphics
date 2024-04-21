library(testthat)
library(MultiModalGraphics)

test_that("ClearScatterplot is properly instantiated", {
  data <- get_clear_scatterplot_df()
  sp <-
    ClearScatterplot(
      data = data,
      timePointColumn = "timePoint",
      logFoldChange = "log2fc",
      timePointLevels = c("T10R1", "T5R1")
    )
  expect_s4_class(sp, "ClearScatterplot")
})

test_that("createPlot modifies the plot object correctly", {
  data <- get_clear_scatterplot_df()
  sp <-
    ClearScatterplot(
      data = data,
      timePointColumn = "timePoint",
      logFoldChange = "log2fc",
      timePointLevels = c("T10R1", "T5R1")
    )
  plot <-
    createPlot(
      sp,
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
  expect_true(!is.null(plot@plot), "Plot should not be NULL after creation")
})

test_that("Errors are thrown for bad input", {
  expect_error(ClearScatterplot(NULL), "Data must be a data frame")
})
