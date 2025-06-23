# tests/testthat/test_multifeaturegrid_plot_basics.R
library(testthat)
library(MultiModalGraphics)
library(ggplot2)

# build a full minimal data.frame
df_full <- data.frame(
  tissue        = factor(rep(c("T1","T2"), each=3)),
  signaling     = factor(rep(c("S1","S2","S3"), 2)),
  Activation_z_score = rnorm(6),
  p             = runif(6, min=0.001, max=0.05),
  number_of_genes   = sample(1:10, 6, replace=TRUE),
  timePoint     = rep(c("Early","Late"), each=3),
  stringsAsFactors = FALSE
)

mg <- CompositeFeatureHeatmap(df_full)

test_that("plot_heatmap returns a ggplot and prints invisibly", {
  # Should return a ggplot object without error
  p <- plot_heatmap(mg)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap errors if pValueColumn contains zeros or negatives", {
  df2 <- df_full
  df2$p[1] <- 0
  mg2 <- CompositeFeatureHeatmap(df2)
  expect_error(
    plot_heatmap(mg2, pValueColumn = "p"),
    "must be strictly > 0"
  )
})

test_that("plot_heatmap errors if asked for non-existent pValueColumn", {
  expect_error(
    plot_heatmap(mg, pValueColumn = "no_such_col"),
    "not found in data"
  )
})

test_that("plot_heatmap errors if columnForNumber missing or negative", {
  df3 <- df_full
  df3$number_of_genes[1] <- -1
  mg3 <- CompositeFeatureHeatmap(df3)
  expect_error(
    plot_heatmap(mg3, columnForNumber = "number_of_genes"),
    "must be ≥ 0"
  )
  expect_error(
    plot_heatmap(mg, columnForNumber = "no_such_col"),
    "not found in data"
  )
})
