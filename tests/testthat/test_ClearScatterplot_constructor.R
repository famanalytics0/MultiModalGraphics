# tests/testthat/test_ThresholdedScatterplot_constructor.R
library(testthat)
library(MultiModalGraphics)

test_that("ThresholdedScatterplot constructor works with valid data", {
  df <- data.frame(
    log2fc     = c(2, -3, 0.5),
    negLog10p  = c(2.0, 2.5, 0.8),
    regulation = c("up", "down", "neutral"),
    SampleType = factor(c("A","B","A")),
    stringsAsFactors = FALSE
  )
  cs <- ThresholdedScatterplot(df)
  expect_s4_class(cs, "ThresholdedScatterplot")
  expect_true(is.data.frame(cs@data))
  expect_identical(levels(cs@data$category), c("down","neutral","up"))
})

test_that("ThresholdedScatterplot errors when required columns are missing or wrong type", {
  df1 <- data.frame( log2fc=1:3, negLog10p=1:3 )
  expect_error(ThresholdedScatterplot(df1), "Missing columns")
  expect_error(ThresholdedScatterplot("not a data.frame"), "is.data.frame")
})

test_that("ThresholdedScatterplot drops NA rows with a warning", {
  df <- data.frame(
    log2fc     = c(1, NA, 2),
    negLog10p  = c(1, 2,  NA),
    regulation = c("up","down","up"),
    SampleType = c("A","A","B"),
    stringsAsFactors = FALSE
  )
  expect_warning(cs <- ThresholdedScatterplot(df), "row\\(s\\) removed")
  expect_equal(nrow(cs@data), 1)
})
