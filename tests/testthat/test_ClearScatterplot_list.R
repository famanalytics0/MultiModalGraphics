library(testthat)
library(MultiModalGraphics)

test_that("list‐of‐matrices constructor single vs multi", {
  # single
  mat1 <- matrix(rnorm(12), nrow = 3)
  colnames(mat1) <- paste0("S",1:4)
  meta1 <- data.frame(Group="G1", SampleType=paste0("S",1:4), row.names=colnames(mat1))
  cs1 <- ThresholdedScatterplot(data = list(one = mat1), meta = list(one = meta1))
  expect_s4_class(cs1, "ThresholdedScatterplot")
  # multi
  mat2 <- matrix(rnorm(12), nrow = 3)
  colnames(mat2) <- paste0("T",1:4)
  meta2 <- data.frame(Group="G1", SampleType=paste0("T",1:4), row.names=colnames(mat2))
  cs2 <- ThresholdedScatterplot(
    data = list(M1 = mat1, M2 = mat2),
    meta = list(M1 = meta1, M2 = meta2)
  )
  expect_s4_class(cs2, "ThresholdedScatterplot")
  # faceting by SampleType appears
  expect_true("timePoint ~ SampleType" %in% as.character(cs2@plot$layers[[1]]$mapping))
})
