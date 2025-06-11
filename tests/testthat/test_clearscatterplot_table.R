library(testthat)
library(MultiModalGraphics)

test_that("matrix+meta constructor works", {
  mat <- matrix(rnorm(20), nrow = 5)
  colnames(mat) <- paste0("S", 1:4)
  meta <- data.frame(
    Group      = rep(c("G1","G2"), each = 2),
    SampleType = paste0("S", 1:4),
    stringsAsFactors = FALSE,
    row.names = paste0("S", 1:4)
  )
  cs <- ClearScatterplot_table(
    expr        = mat,
    meta        = meta,
    groupColumn = "Group",
    sampleType  = "SampleType"
  )
  expect_s4_class(cs, "ClearScatterplot")
  expect_equal(sort(unique(cs@data$SampleType)), sort(unique(meta$SampleType)))
})
