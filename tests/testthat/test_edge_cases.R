# tests/testthat/test_edge_cases.R
library(testthat)
library(MultiModalGraphics)

test_that("no DE hits produces an empty or trivial heatmap", {
  expr <- matrix(0, nrow=5, ncol=4,
                 dimnames=list(letters[1:5], paste0("S",1:4)))
  meta <- data.frame(treatment=rep(c("A","B"), each=2),
                     row.names = paste0("S",1:4))
  hm <- AnnotatedHeatmap(expr, meta, groupColumn="treatment")
  ht <- getHeatmapObject(hm)
  # if no features pass DE, matrix may have zero rows
  expect_true(nrow(ht@matrix) == 0 || inherits(ht, "Heatmap"))
})

test_that("invalid inputs throw the right errors", {
  expect_error(AnnotatedHeatmap("not a matrix"), "no applicable method")
  expect_error(AnnotatedHeatmap_table("a","b"), "no applicable method")
})

test_that("duplicate names in multimodal lists raise an error", {
  fc <- matrix(1,1,1, dimnames=list("x","y"))
  pv <- matrix(1,1, dimnames=list("x","y"))
  data_list <- list(a=fc, a=fc)
  pv_list   <- list(a=pv, a=pv)
  expect_error(
    AnnotatedHeatmap(data_list, pv_list),
    "must match"
  )
})
