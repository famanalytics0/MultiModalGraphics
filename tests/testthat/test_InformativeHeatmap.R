library(testthat)
library(ComplexHeatmap)
library(MultiModalGraphics)

test_that("AnnotatedHeatmap initializes correctly", {
  data <- matrix(rnorm(100), ncol = 10)
  heatmap <- AnnotatedHeatmap(data)
  expect_s4_class((heatmap), "AnnotatedHeatmap")
})

test_that("AnnotatedHeatmap handles invalid input", {
  data <- "not a matrix"
  expect_error(AnnotatedHeatmap(data))
})

test_that("updateLayerFun updates layer function correctly", {
  data <- matrix(rnorm(100), ncol = 10)
  heatmap <- AnnotatedHeatmap(data)
  new_layer_fun <- function(j, i, x, y, w, h, fill) {
    grid::grid.points(x, y, pch = 21, size = unit(2, "mm"), gp = grid::gpar(col = "red"))
  }
  updated_heatmap <- updateLayerFun(heatmap, new_layer_fun)
  # Inspect if the new function is actually being used
  # This might require your function to allow introspection or store metadata about the last applied function
  expect_equal(updated_heatmap@params$layer_fun, new_layer_fun)
})

test_that("getHeatmapObject returns a Heatmap object", {
  data <- matrix(rnorm(100), ncol = 10)
  heatmap <- AnnotatedHeatmap(data)
  heatmap_obj <- getHeatmapObject(heatmap)
  expect_s4_class(heatmap_obj, "Heatmap")
})

