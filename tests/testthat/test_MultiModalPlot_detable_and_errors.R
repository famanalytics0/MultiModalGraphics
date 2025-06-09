# tests/testthat/test_MultiModalPlot_detable_and_errors.R
library(testthat)
library(MultiModalGraphics)

test_that("DE‐table → volcano only, and errors on heatmap", {
  df <- data.frame(
    log2fc    = rnorm(10),
    negLog10p = runif(10,1,5),
    regulation= sample(c("up","down"),10,TRUE),
    SampleType= sample(c("A","B"),10,TRUE),
    stringsAsFactors=FALSE
  )
  p <- MultiModalPlot(
    inputs       = list(M2=df),
    panel_type   = "volcano"
  )
  expect_s3_class(p, "ggplot")
  expect_error(
    MultiModalPlot(inputs=list(M2=df), panel_type="heatmap"),
    "DE-tables only support panel_type = 'volcano'"
  )
})

test_that("unrecognized inputs throw", {
  expect_error(
    MultiModalPlot(inputs = list(bad=42)),
    "Unrecognized input type for 'bad'"
  )
})
