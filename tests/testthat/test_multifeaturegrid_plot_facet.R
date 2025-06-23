# tests/testthat/test_multifeaturegrid_plot_facet.R
library(testthat)
library(MultiModalGraphics)
library(ggplot2)

# timePoint all NA → skip facet with a warning
df_na_facet <- df_full
df_na_facet$timePoint <- NA
mg_na <- CompositeFeatureHeatmap(df_na_facet)

test_that("plot_heatmap skips facet when facet column is all NA", {
  expect_warning(
    p_na <- plot_heatmap(mg_na, independantVariable = "timePoint"),
    "facet will be skipped"
  )
  expect_s3_class(p_na, "ggplot")
})

test_that("plot_heatmap without specifying facet still works", {
  p_none <- plot_heatmap(mg, independantVariable = "nonexistent_col")
  expect_s3_class(p_none, "ggplot")
})
