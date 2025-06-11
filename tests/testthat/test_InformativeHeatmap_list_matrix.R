library(testthat)
library(MultiModalGraphics)

test_that("InformativeHeatmap(list of matrices) single vs multi", {
  set.seed(2)
  # single
  mat1 <- matrix(rnorm(12), nrow=3)
  rownames(mat1) <- paste0("g",1:3)
  colnames(mat1) <- paste0("s",1:4)
  meta1 <- data.frame(
    Group = rep(1:2, each=2),
    row.names = colnames(mat1),
    stringsAsFactors = FALSE
  )
  hm1 <- InformativeHeatmap(
    data        = list(M1=mat1),
    meta        = list(M1=meta1),
    groupColumn = "Group"
  )
  expect_s4_class(hm1, "InformativeHeatmap")
  expect_true(hm1@params$single)

  # multi
  mat2 <- matrix(rnorm(12), nrow=3)
  rownames(mat2) <- paste0("h",1:3)
  colnames(mat2) <- paste0("t",1:4)
  meta2 <- data.frame(
    Group = rep(c("X","Y"), each=2),
    row.names = colnames(mat2),
    stringsAsFactors = FALSE
  )
  hm2 <- InformativeHeatmap(
    data        = list(M1=mat1, M2=mat2),
    meta        = list(M1=meta1, M2=meta2),
    groupColumn = "Group"
  )
  expect_s4_class(hm2, "InformativeHeatmap")
  expect_false(isTRUE(hm2@params$single))
})
