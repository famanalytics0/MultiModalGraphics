library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)

test_that("AnnotatedHeatmap(matrix, meta) works for single modality", {
  # toy expression matrix
  set.seed(1)
  mat <- matrix(rnorm(20), nrow=5)
  rownames(mat) <- paste0("g", 1:5)
  colnames(mat) <- paste0("s", 1:4)
  meta <- data.frame(
    Group = rep(c("A","B"), each=2),
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )
  hm <- AnnotatedHeatmap(
    data        = mat,
    meta        = meta,
    groupColumn = "Group",
    heatmap_data_scale = "expression"
  )
  expect_s4_class(hm, "AnnotatedHeatmap")
  # heatmap object exists
  ht_obj <- getHeatmapObject(hm)
  expect_true(inherits(ht_obj, "Heatmap"))
  # params recorded
  expect_equal(hm@params$groupColumn, "Group")
  expect_equal(hm@params$heatmap_data_scale, "expression")
})
