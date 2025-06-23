# tests/testthat/test_table_constructor.R
library(testthat)
library(MultiModalGraphics)

test_that("table-based constructor returns AnnotatedHeatmap", {
  fc <- matrix(rnorm(20), nrow=5, dimnames=list(letters[1:5], paste0("G",1:4)))
  pv <- matrix(runif(20), nrow=5, dimnames=list(letters[1:5], paste0("G",1:4)))
  hm <- AnnotatedHeatmap_table(fc, pv)
  expect_s4_class(hm, "AnnotatedHeatmap")
  expect_s4_class(getHeatmapObject(hm), "Heatmap")
})

test_that("table-based constructor errors on mismatched dims", {
  fc <- matrix(rnorm(6), nrow=3)
  pv <- matrix(runif(8), nrow=4)
  expect_error(
    AnnotatedHeatmap_table(fc, pv),
    "non-conformable arrays|dimensions"
  )
})
