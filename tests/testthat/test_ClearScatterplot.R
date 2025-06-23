library(testthat)
library(MultiModalGraphics)

test_that("Plain DE‐data constructor works", {
  df <- data.frame(
    log2fc     = c(-2, 0, 2),
    negLog10p  = c(3, 1, 4),
    regulation = c("down","neutral","up"),
    SampleType = rep("A", 3),
    stringsAsFactors = FALSE
  )
  cs <- ThresholdedScatterplot(df)
  expect_s4_class(cs, "ThresholdedScatterplot")
  expect_equal(levels(cs@data$category), c("down","neutral","up"))
  # make sure plotting slot gets filled
  cs <- createPlot(cs)
  expect_true(inherits(cs@plot, "ggplot"))
})
