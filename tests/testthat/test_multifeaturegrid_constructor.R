# tests/testthat/test_multifeaturegrid_constructor.R
library(testthat)
library(MultiModalGraphics)  # replace with your package

test_that("constructor succeeds with minimal required columns", {
  df <- data.frame(
    tissue = factor(c("T1","T2")),
    signaling = factor(c("S1","S2")),
    Activation_z_score = c(-1.5, 2.3),
    stringsAsFactors = FALSE
  )
  mg <- CompositeFeatureHeatmap(df)
  expect_s4_class(mg, "CompositeFeatureHeatmap")
  expect_identical(slot(mg, "breaks"), seq(-1,1,0.5))
})

test_that("constructor errors when required columns are missing", {
  df1 <- data.frame(
    tissue = factor(c("T1","T2")),
    Activation_z_score = c(0,1)
  )
  expect_error(
    CompositeFeatureHeatmap(df1),
    "`data` is missing required column\\(s\\): signaling"
  )

  df2 <- data.frame(
    tissue = factor(c("T1","T2")),
    signaling = factor(c("S1","S2")),
    Activation_z_score = as.character(c("a","b")),
    stringsAsFactors = FALSE
  )
  expect_error(
    CompositeFeatureHeatmap(df2),
    "Column `Activation_z_score` must be numeric"
  )
})

test_that("constructor enforces strictly increasing breaks", {
  df <- data.frame(
    tissue = factor(c("T1","T2")),
    signaling = factor(c("S1","S2")),
    Activation_z_score = c(0,0),
    stringsAsFactors = FALSE
  )
  expect_error(
    CompositeFeatureHeatmap(df, breaks = c(0,0,1)),
    "`breaks` must be strictly increasing"
  )
  expect_error(
    CompositeFeatureHeatmap(df, breaks = "not numeric"),
    "`breaks` must be a numeric vector"
  )
})

test_that("constructor errors on invalid palette", {
  df <- data.frame(
    tissue = factor(c("T1","T2")),
    signaling = factor(c("S1","S2")),
    Activation_z_score = c(0,0),
    stringsAsFactors = FALSE
  )
  expect_error(
    CompositeFeatureHeatmap(df, color_palette = "NotAPalette"),
    "`color_palette` must be a valid RColorBrewer palette name"
  )
})

