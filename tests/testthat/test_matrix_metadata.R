# tests/testthat/test_matrix_metadata.R
library(testthat)
library(MultiModalGraphics) # or library(AnnotatedHeatmap)

test_that("matrix+metadata constructor errors when groupColumn is missing", {
  expr <- matrix(rnorm(20), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:4)))
  meta <- data.frame(treatment = rep(c("A","B"), each=2),
                     row.names = colnames(expr))
  expect_error(
    AnnotatedHeatmap(expr, meta),
    "must supply `groupColumn`"
  )
})

test_that("matrix+metadata constructor succeeds with valid groupColumn", {
  expr <- matrix(rnorm(20), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:4)))
  meta <- data.frame(treatment = rep(c("A","B"), each=2),
                     row.names = colnames(expr))
  hm <- AnnotatedHeatmap(expr, meta, groupColumn = "treatment")
  expect_s4_class(hm, "AnnotatedHeatmap")
  expect_s4_class(getHeatmapObject(hm), "Heatmap")
})
